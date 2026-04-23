import json
import os
import sys
import time
import copy
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Union, Optional
import warnings 
import argparse

import matplotlib
import numpy as np
import pandas as pd
import torch
import wandb
from scipy.sparse import csr_matrix
from torch import nn
from torch.nn import functional as F
from torchtext.vocab import Vocab
from torchtext._torchtext import Vocab as VocabPybind
from torch_geometric.loader import DataLoader
from torch.utils.data import ConcatDataset
from gears import PertData, GEARS
from gears.inference import compute_metrics, deeper_analysis, non_dropout_analysis
from gears.utils import create_cell_graph_dataset_for_prediction

PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
sys.path.insert(0, os.path.join(PATH, "scripts/scGPT"))

import scgpt as scg
from scgpt.model import TransformerGenerator
from scgpt.loss import masked_mse_loss, criterion_neg_log_bernoulli, masked_relative_error
from scgpt.tokenizer import tokenize_batch, pad_batch, tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed, map_raw_id_to_vocab_id, compute_perturbation_metrics


def train(model: nn.Module, train_loader: torch.utils.data.DataLoader, device: torch.device, configs: dict, logger) -> None:
    """
    Train the model for one epoch.
    """
    
    required_keys = [
        "include_zero_gene",
        "max_seq_len",
        "gene_ids",
        "amp",
        "CLS",
        "CCE",
        "MVC",
        "ECS",
        "log_interval",
        "criterion",
        "optimizer",
        "scaler",
        "scheduler",
        "epoch",
        "pert_list",
        "gene_map"
    ]

    for key in required_keys:
        assert key in configs, f"Missing required hyperparameter: '{key}'"
        
    include_zero_gene = configs["include_zero_gene"]
    max_seq_len = configs["max_seq_len"]
    gene_ids = configs["gene_ids"]
    amp = configs["amp"]
    CLS = configs["CLS"]
    CCE = configs["CCE"]
    MVC = configs["MVC"]
    ECS = configs["ECS"]
    log_interval = configs["log_interval"]
    scheduler = configs["scheduler"]
    criterion = configs["criterion"]
    optimizer = configs["optimizer"]
    scaler = configs["scaler"]
    epoch = configs["epoch"]
    pert_names = configs["pert_list"]  
    gene_map = configs["gene_map"] 
    
    model.train()
    total_loss, total_mse = 0.0, 0.0
    start_time = time.time()

    num_batches = len(train_loader)
    for batch, batch_data in enumerate(train_loader):
        batch_size, n_genes = batch_data.y.shape
        batch_data.to(device)
        x: torch.Tensor = batch_data.x  # (batch_size * n_genes, 1)
        ori_gene_values = x[:, 0].view(batch_size, n_genes)

        # Track which cells to keep (scGPT can only predict perturbations for which gene is also sequenced)
        valid_cells = []
        # Adding perturbation vector
        # https://github.com/snap-stanford/GEARS/blob/c7ca19cbcd6f4da3030d0ebc90b2c2cd0b47a8d8/gears/pertdata.py#L360-L364
        pert_feats = torch.zeros(x.shape, device=device)
        pert_idx = batch_data.pert_idx
        if pert_idx is not None:
            for i, p in enumerate(pert_idx):
                perturbed_gene = pert_names[int(np.abs(p))] # maps pert_idx to gene_name
                # Check if gene is in gene_map
                if perturbed_gene not in gene_map:
                    continue
                gene_position_in_adata = gene_map[perturbed_gene]  # maps to adata.var.gene_names
                pert_feats[n_genes*i + gene_position_in_adata, 0] = 1
                valid_cells.append(i)
        # If no valid cells, skip this batch
        if len(valid_cells) == 0:
            logger.warning("No valid cells in batch. Skipping batch.")
            continue
            # Filter to keep only valid cells
        if len(valid_cells) < batch_size:
            # Create mask for valid cells
            valid_mask = torch.zeros(batch_size, dtype=torch.bool, device=device)
            valid_mask[valid_cells] = True

            # Reshape to (batch_size, n_genes) for easier filtering
            x_reshaped = x.view(batch_size, n_genes, -1)  # (batch_size, n_genes, 1)
            pert_feats_reshaped = pert_feats.reshape(batch_size, n_genes, -1)

            # Filter valid cells
            x_reshaped = x_reshaped[valid_mask]  # (valid_batch_size, n_genes, 1)
            pert_feats_reshaped = pert_feats_reshaped[valid_mask]  # (valid_batch_size, n_genes, 1)
            batch_data.y = batch_data.y[valid_mask]  # (valid_batch_size, n_genes)

            # Flatten back
            valid_batch_size = len(valid_cells)
            x = x_reshaped.reshape(valid_batch_size * n_genes, -1)
            pert_feats = pert_feats_reshaped.reshape(valid_batch_size * n_genes, -1)
            ori_gene_values = x[:, 0].view(valid_batch_size, n_genes)
            batch_size = valid_batch_size
            
        x = torch.cat((x, pert_feats), dim = 1).to(torch.float32)
        pert_flags = x[:, 1].long().view(batch_size, n_genes)
        target_gene_values = batch_data.y  # (batch_size, n_genes)

        if include_zero_gene in ["all", "batch-wise"]:
            if include_zero_gene == "all":
                input_gene_ids = torch.arange(n_genes, device=device, dtype=torch.long)
            else:
                input_gene_ids = (
                    ori_gene_values.nonzero()[:, 1].flatten().unique().sort()[0]
                )
            # sample input_gene_id
            if len(input_gene_ids) > max_seq_len:
                input_gene_ids = torch.randperm(len(input_gene_ids), device=device)[
                    :max_seq_len
                ]
            input_values = ori_gene_values[:, input_gene_ids]
            input_pert_flags = pert_flags[:, input_gene_ids]
            target_values = target_gene_values[:, input_gene_ids]

            mapped_input_gene_ids = map_raw_id_to_vocab_id(input_gene_ids, gene_ids)
            mapped_input_gene_ids = mapped_input_gene_ids.repeat(batch_size, 1)

            # src_key_padding_mask = mapped_input_gene_ids.eq(vocab[pad_token])
            src_key_padding_mask = torch.zeros_like(
                input_values, dtype=torch.bool, device=device
            )

        with torch.cuda.amp.autocast(enabled=amp):
            output_dict = model(
                mapped_input_gene_ids,
                input_values,
                input_pert_flags,
                src_key_padding_mask=src_key_padding_mask,
                CLS=CLS,
                CCE=CCE,
                MVC=MVC,
                ECS=ECS,
            )
            output_values = output_dict["mlm_output"]

            masked_positions = torch.ones_like(
                input_values, dtype=torch.bool
            )  # Use all
            loss = loss_mse = criterion(output_values, target_values, masked_positions)

        model.zero_grad()
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("always")
            torch.nn.utils.clip_grad_norm_(
                model.parameters(),
                1.0,
                error_if_nonfinite=False if scaler.is_enabled() else True,
            )
            if len(w) > 0:
                logger.warning(
                    f"Found infinite gradient. This may be caused by the gradient "
                    f"scaler. The current scale is {scaler.get_scale()}. This warning "
                    "can be ignored if no longer occurs after autoscaling of the scaler."
                )
        scaler.step(optimizer)
        scaler.update()
        
        wandb.log({"train/loss": loss.item()})
        # torch.cuda.empty_cache()

        total_loss += loss.item()
        total_mse += loss_mse.item()
        if batch % log_interval == 0 and batch > 0:
            lr = scheduler.get_last_lr()[0]
            ms_per_batch = (time.time() - start_time) * 1000 / log_interval
            cur_loss = total_loss / log_interval
            cur_mse = total_mse / log_interval
            # ppl = math.exp(cur_loss)
            logger.info(
                f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
                f"lr {lr:05.4f} | ms/batch {ms_per_batch:5.2f} | "
                f"loss {cur_loss:5.2f} | mse {cur_mse:5.2f} |"
            )
            total_loss = 0
            total_mse = 0
            start_time = time.time()


def eval_perturb(
    loader: DataLoader, model: TransformerGenerator, device: torch.device, configs: dict, logger
) -> Dict:
    """
    Run model in inference mode using a given data loader
    """

    model.eval()
    model.to(device)
    pert_cat = []
    pred = []
    truth = []
    pred_de = []
    truth_de = []
    results = {}
    logvar = []
    
    required_keys = [
        "include_zero_gene",
        "gene_ids",
        "pert_list",
        "gene_map"
    ]

    for key in required_keys:
        assert key in configs, f"Missing required hyperparameter: '{key}'"
                                                  
    include_zero_gene = configs["include_zero_gene"]
    gene_ids = configs["gene_ids"]
    pert_names = configs["pert_list"]  
    gene_map = configs["gene_map"] 
                                                  
    for itr, batch in enumerate(loader):
        batch.to(device)
                                                  
        x: torch.Tensor = batch.x  # (batch_size * n_genes, 1)
        batch_size, n_genes = batch.y.shape
        # Track which cells to keep (scGPT can only predict perturbations for which gene is also sequenced)
        valid_cells = []
        # Adding perturbation vector
        # https://github.com/snap-stanford/GEARS/blob/c7ca19cbcd6f4da3030d0ebc90b2c2cd0b47a8d8/gears/pertdata.py#L360-L364
        pert_feats = torch.zeros(x.shape, device=device)
        pert_idx = batch.pert_idx[0]
        if pert_idx is not None:
            for i, p in enumerate(pert_idx):
                perturbed_gene = pert_names[int(np.abs(p))] # maps pert_idx to gene_name
                if perturbed_gene not in gene_map:
                    continue
                gene_position_in_adata = gene_map[perturbed_gene]  # maps to adata.var.gene_names
                pert_feats[n_genes*i + gene_position_in_adata, 0] = 1
                valid_cells.append(i)
        # If no valid cells, skip this batch
        if len(valid_cells) == 0:
            logger.warning(f"No valid cells in batch {itr}. Skipping batch.")
            continue 
        if len(valid_cells) < batch_size:
            # Create mask for valid cells
            valid_mask = torch.zeros(batch_size, dtype=torch.bool, device=device)
            valid_mask[valid_cells] = True

            # Reshape to (batch_size, n_genes) for easier filtering
            x_reshaped = x.view(batch_size, n_genes, -1)  # (batch_size, n_genes, 1)
            pert_feats_reshaped = pert_feats.reshape(batch_size, n_genes, -1)

            # Filter valid cells
            x_reshaped = x_reshaped[valid_mask]  # (valid_batch_size, n_genes, 1)
            pert_feats_reshaped = pert_feats_reshaped[valid_mask]  # (valid_batch_size, n_genes, 1)
            batch.y = batch.y[valid_mask]  # (valid_batch_size, n_genes)
            batch.pert_idx = [batch.pert_idx[i] for i in valid_cells]
            batch.de_idx = [batch.de_idx[i] for i in valid_cells]
            batch.pert = [batch.pert[i] for i in valid_cells]
            
            # Flatten back
            valid_batch_size = len(valid_cells)
            x = x_reshaped.reshape(valid_batch_size * n_genes, -1)
            pert_feats = pert_feats_reshaped.reshape(valid_batch_size * n_genes, -1)
            batch_size = valid_batch_size
        
        pert_cat.extend(batch.pert)
        pert_feats = torch.tensor(pert_feats, device=device)
        batch.x = torch.cat((x, pert_feats), dim = 1).to(torch.float32)
                                                  
        with torch.no_grad():
            p = model.pred_perturb(
                batch,
                include_zero_gene=include_zero_gene,
                gene_ids=gene_ids,
            )
            t = batch.y
            pred.extend(p.cpu())
            truth.extend(t.cpu())

            # Differentially expressed genes
            for itr, de_idx in enumerate(batch.de_idx):
                pred_de.append(p[itr, de_idx])
                truth_de.append(t[itr, de_idx])

    # all genes
    results["pert_cat"] = np.array(pert_cat)
    pred = torch.stack(pred)
    truth = torch.stack(truth)
    results["pred"] = pred.detach().cpu().numpy().astype(np.float)
    results["truth"] = truth.detach().cpu().numpy().astype(np.float)

    pred_de = torch.stack(pred_de)
    truth_de = torch.stack(truth_de)
    results["pred_de"] = pred_de.detach().cpu().numpy().astype(np.float)
    results["truth_de"] = truth_de.detach().cpu().numpy().astype(np.float)

    return results

def main():
    parser = argparse.ArgumentParser(description="A script that fine-tunes a scGPT model")

    parser.add_argument(
        '--data_name',
        type=str,
        required=True,
        help='The name of the data to be processed.'
    )
    parser.add_argument(
        '--split_type',
        type=str,
        required=True,
        help='The type of split: 1/10th or 9/10th'
    )
    parser.add_argument(
        '--split_num',
        type=int,
        required=True,
        help='The specific split number of the data to be processed (e.g., for cross-validation).'
    )
    
    args = parser.parse_args()
    print(f"Processing {args.data_name} on split {args.split_num} with {args.split_type} split type.")

    hyperparameter_defaults = dict(
        seed = 0,
        split = args.split_num,
        split_type = args.split_type,
        data_name = args.data_name,
        load_model = os.path.join(PATH, "scripts/scGPT/save/scGPT_human"),
        max_seq_len = 1536,
        lr = 1e-4,  # or 1e-4
        batch_size = 64,
        eval_batch_size = 64,
        epochs = 15,
        schedule_interval = 1,
        early_stop = 10,
        embsize = 512,  # embedding dimension
        d_hid = 512,  # dimension of the feedforward network model in nn.TransformerEncoder
        nlayers = 12,  # number of nn.TransformerEncoderLayer in nn.TransformerEncoder
        nhead = 8,  # number of heads in nn.MultiheadAttention
        n_layers_cls = 3,
        dropout = 0,  # dropout probability
        log_interval = 100,
        use_fast_transformer = True,  # whether to use fast transformer
        MLM = True,  # whether to use masked language modeling, currently it is always on.
        CLS = False,  # celltype classification objective
        CCE = False,  # Contrastive cell embedding objective
        MVC = False,  # Masked value prediction for cell embedding
        ECS = False,  # Elastic cell similarity objective
        amp = True,
        include_zero_gene = "all"
    )

    run = wandb.init(
        project="sigrecon",
        name=f"scgpt_{args.data_name}_{args.split_type}_{args.split_num}",
        config=hyperparameter_defaults,
        reinit=True,
        settings=wandb.Settings(start_method="fork"),
    )

    config = wandb.config
    set_seed(config.seed)

    # settings for data prcocessing
    pad_token = "<pad>"
    special_tokens = [pad_token, "<cls>", "<eoc>"]
    pad_value = 0  # for padding values
    pert_pad_id = 0
    max_seq_len = config.max_seq_len

    load_param_prefixs = [
        "encoder",
        "value_encoder",
        "transformer_encoder",
    ]
    
    # settings for optimizer
    batch_size = config.batch_size
    eval_batch_size = config.batch_size
    schedule_interval = config.schedule_interval

    # settings for the model
    embsize = config.embsize  # embedding dimension
    d_hid = config.embsize  # dimension of the feedforward network in transformer
    nlayers = config.nlayers  # number of transformer layers
    nhead = config.nhead  # number of heads in nn.MultiheadAttention
    n_layers_cls = 3

    # dataset and evaluation choices
    split = "no_test"

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    save_dir = Path(f"./save/dev_perturb_{args.data_name}-{args.split_num}-{time.strftime('%b%d-%H-%M')}/")
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"saving to {save_dir}")

    logger = scg.logger
    scg.utils.add_file_handler(logger, save_dir / "run.log")
    # log running date and current git commit
    logger.info(f"Running on {time.strftime('%Y-%m-%d %H:%M:%S')}")    

    # Loading Pert Data
    pert_data_source = PertData(os.path.join(PATH, "data/scgpt"))
    pert_data_target = PertData(os.path.join(PATH, "data/scgpt"))
    if args.data_name == 'k562':
        # pert_data_source.load(data_path = os.path.join(PATH, "data/scgpt/k562"))
        pert_data_source.load(data_path = os.path.join(PATH, "data/scgpt/k562_shared"))
        pert_data_target.load(data_path = os.path.join(PATH, "data/scgpt/rpe1_shared"))
    elif args.data_name == "rpe1":
        # pert_data_source.load(data_path = os.path.join(PATH, "data/scgpt/rpe1"))
        pert_data_source.load(data_path = os.path.join(PATH, "data/scgpt/rpe1_shared"))
        pert_data_target.load(data_path = os.path.join(PATH, "data/scgpt/k562_shared"))
        
    pert_data_source.prepare_split(split = "no_split", seed = 42) 
    pert_data_source_dataloader = pert_data_source.get_dataloader(batch_size=batch_size, test_batch_size=eval_batch_size)
    pert_data_target.prepare_split(split = "no_split", seed = 42) 
    pert_data_target_dataloader=pert_data_target.get_dataloader(batch_size=batch_size, test_batch_size=batch_size)
    
    # Concatenate the datasets (target control + 1/10 or 9/10 of target unique perturbations)
    splits_df = pd.read_csv(os.path.join(PATH, "data/sigs/perturb-seq/pb_splits.csv")) # splits df is all in rpe1_shared and k562_shared
    split_col = f"split_{args.split_num}"
    
    # Get unique perturbations from the target dataset
    unique_perturbations = np.unique([data["pert_idx"][0] for data in pert_data_target_dataloader["test_loader"].dataset])

    # Extract perturbation names and idx from target
    pert_idx_to_name = {data["pert_idx"][0]: data["pert"].replace("+ctrl","") for data in pert_data_target_dataloader["test_loader"].dataset}
    pert_name_to_idx = {v: k for k, v in pert_idx_to_name.items()}

    # Get the split column name
    split_col = f"split_{args.split_num}"

    # Get perturbations in the eval split (TRUE values = held out 1/10th)
    if args.split_type == "90th":
        train_pb_names = splits_df[splits_df[split_col] == False]["pb"].values
        eval_pb_names = splits_df[splits_df[split_col] == True]["pb"].values
    elif args.split_type == "10th":
        train_pb_names = splits_df[splits_df[split_col] == True]["pb"].values
        eval_pb_names = splits_df[splits_df[split_col] == False]["pb"].values

    print(f"Using split {args.split_num}: {len(eval_pb_names)} eval perturbations, {len(train_pb_names)} train perturbations")

    # Convert names to indices (filter to only those present in the dataset)
    eval_perturbations = [pert_name_to_idx[name] for name in eval_pb_names if name in pert_name_to_idx]
    train_perturbations = [pert_name_to_idx[name] for name in train_pb_names if name in pert_name_to_idx]

    print(f"Found {len(eval_perturbations)} eval and {len(train_perturbations)} train perturbations in dataset")

    # Split the target dataset based on predefined splits
    train_data_split = [data for data in pert_data_target_dataloader["test_loader"].dataset 
                        if data["pert_idx"][0] in train_perturbations]
    eval_data_split = [data for data in pert_data_target_dataloader["test_loader"].dataset 
                       if data["pert_idx"][0] in eval_perturbations]

    # Concatenate the datasets (All source, 90% target)
    combined_dataset = ConcatDataset([pert_data_source_dataloader["test_loader"].dataset, train_data_split])

    # Create DataLoaders
    combined_train_dataloader = DataLoader(combined_dataset, batch_size=batch_size, shuffle=True)
    eval_dataloader = DataLoader(eval_data_split, batch_size=batch_size, shuffle=True)
    
    # Adding pert_list and gene_map
    hyperparameter_defaults["pert_list"] = pert_data_source.pert_names # List of perts defined by set_pert_genes, pert_idx indexes into this
    hyperparameter_defaults["gene_map"] = pert_data_source.node_map # Dictionary mapping gene names to gene_var index in adata
    
    if config.load_model is not None:
        model_dir = Path(config.load_model)
        model_config_file = model_dir / "args.json"
        model_file = model_dir / "best_model.pt"
        vocab_file = model_dir / "vocab.json"

        vocab = GeneVocab.from_file(vocab_file)
        for s in special_tokens:
            if s not in vocab:
                vocab.append_token(s)

        pert_data_source.adata.var["id_in_vocab"] = [
            1 if gene in vocab else -1 for gene in pert_data_source.adata.var["gene_name"]
        ]
        gene_ids_in_vocab = np.array(pert_data_source.adata.var["id_in_vocab"])
        logger.info(
            f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
            f"in vocabulary of size {len(vocab)}."
        )
        genes = pert_data_source.adata.var["gene_name"].tolist()

        # model
        with open(model_config_file, "r") as f:
            model_configs = json.load(f)
        logger.info(
            f"Resume model from {model_file}, the model args will override the "
            f"config {model_config_file}."
        )
        embsize = model_configs["embsize"]
        nhead = model_configs["nheads"]
        d_hid = model_configs["d_hid"]
        nlayers = model_configs["nlayers"]
        n_layers_cls = model_configs["n_layers_cls"]
    else:
        genes = pert_data_source.adata.var["gene_name"].tolist()
        vocab = Vocab(
            VocabPybind(genes + special_tokens, None)
        )  # bidirectional lookup [gene <-> int]
    vocab.set_default_index(vocab["<pad>"])

    #gene ids relative to the transformer vocab
    gene_ids = np.array(
        [vocab[gene] if gene in vocab else vocab["<pad>"] for gene in genes], dtype=int
    )
    hyperparameter_defaults["gene_ids"] = gene_ids
    
    # Load Model
    ntokens = len(vocab)  # size of vocabulary
    model = TransformerGenerator(
        ntokens,
        embsize,
        nhead,
        d_hid,
        nlayers,
        nlayers_cls=n_layers_cls,
        n_cls=1,
        vocab=vocab,
        dropout=config.dropout,
        pad_token=pad_token,
        pad_value=pad_value,
        pert_pad_id=pert_pad_id,
        use_fast_transformer=config.use_fast_transformer,
    )
    if load_param_prefixs is not None and config.load_model is not None:
        # only load params that start with the prefix
        model_dict = model.state_dict()
        pretrained_dict = torch.load(model_file)
        pretrained_dict = {
            k: v
            for k, v in pretrained_dict.items()
            if any([k.startswith(prefix) for prefix in load_param_prefixs])
        }
        for k, v in pretrained_dict.items():
            logger.info(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)
    elif config.load_model is not None:
        try:
            model.load_state_dict(torch.load(model_file))
            logger.info(f"Loading all model params from {model_file}")
        except:
            # only load params that are in the model and match the size
            model_dict = model.state_dict()
            pretrained_dict = torch.load(model_file)
            pretrained_dict = {
                k: v
                for k, v in pretrained_dict.items()
                if k in model_dict and v.shape == model_dict[k].shape
            }
            for k, v in pretrained_dict.items():
                logger.info(f"Loading params {k} with shape {v.shape}")
            model_dict.update(pretrained_dict)
            model.load_state_dict(model_dict)
    model.to(device)
    wandb.watch(model)

    hyperparameter_defaults["criterion"] = masked_mse_loss
    hyperparameter_defaults["criterion_cls"] = nn.CrossEntropyLoss()
    hyperparameter_defaults["optimizer"] = torch.optim.Adam(model.parameters(), lr=config.lr)
    hyperparameter_defaults["scheduler"] = torch.optim.lr_scheduler.StepLR(hyperparameter_defaults["optimizer"], schedule_interval, gamma=0.9)
    hyperparameter_defaults["scaler"] = torch.cuda.amp.GradScaler(enabled=config.amp)
    
    # Training Loop
    best_val_loss = float("inf")
    best_val_corr = 0
    best_model = None
    patience = 0

    for epoch in range(1, config.epochs + 1):
        epoch_start_time = time.time()
        # Train is 90% of data, valid is 10% of data
        # train_loader = pert_data.dataloader["train_loader"]
        train_loader = combined_train_dataloader
        # valid_loader = pert_data_target.dataloader["val_loader"]
        valid_loader = eval_dataloader
        
        hyperparameter_defaults["epoch"] = epoch
        
        train(
            model,
            train_loader,
            device,
            hyperparameter_defaults,
            logger
        )

        val_res = eval_perturb(valid_loader, 
                               model, 
                               device, 
                               hyperparameter_defaults,
                               logger)
        val_metrics = compute_perturbation_metrics(
            val_res, pert_data_target.adata[pert_data_target.adata.obs["condition"] == "ctrl"]
        )
        logger.info(f"val_metrics at epoch {epoch}: ")
        logger.info(val_metrics)
        wandb.log({f"valid/{k}": v for k, v in val_metrics.items()})

        elapsed = time.time() - epoch_start_time
        logger.info(f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | ")

        val_score = val_metrics["pearson"]
        if val_score > best_val_corr:
            best_val_corr = val_score
            best_model = copy.deepcopy(model)
            logger.info(f"Best model with score {val_score:5.4f}")
            patience = 0
        else:
            patience += 1
            if patience >= config.early_stop:
                logger.info(f"Early stop at epoch {epoch}")
                break

        torch.save(
            model.state_dict(),
            save_dir / f"model_{epoch}.pt",
        )

        hyperparameter_defaults["scheduler"].step()

    # Save Best Model    
    torch.save(best_model.state_dict(), save_dir / "best_model.pt")

    wandb.finish()

if __name__ == "__main__":
    main()