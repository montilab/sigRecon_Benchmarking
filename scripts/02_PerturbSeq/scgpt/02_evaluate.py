import json
import os
import sys
import time
import argparse
import glob
from pathlib import Path
from typing import Dict
import numpy as np
import torch
import warnings
import anndata
import pandas as pd

# Assuming MLAB environment variable is set for the base path
PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
sys.path.insert(0, os.path.join(PATH, "scripts/scGPT"))

import scgpt as scg
from scgpt.model import TransformerGenerator
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from gears import PertData
from torch_geometric.loader import DataLoader
from torch.utils.data import ConcatDataset

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
        "gene_names",
        "pert_list",
        "gene_map"
    ]

    for key in required_keys:
        assert key in configs, f"Missing required hyperparameter: '{key}'"
                                                  
    include_zero_gene = configs["include_zero_gene"]
    gene_ids = configs["gene_ids"]
    gene_names = configs["gene_names"]
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

            # # Differentially expressed genes
            # for itr, de_idx in enumerate(batch.de_idx):
            #     pred_de.append(p[itr, de_idx])
            #     truth_de.append(t[itr, de_idx])

    # all genes
    results["pert_cat"] = np.array(pert_cat)
    pred = torch.stack(pred)
    # truth = torch.stack(truth)
    results["pred"] = pred.detach().cpu().numpy().astype(np.float32)
    # results["truth"] = truth.detach().cpu().numpy().astype(np.float32)
    
    # Create and return AnnData object
    obs_df = pd.DataFrame({
        'perturbation': results['pert_cat'].tolist()
    })
    var_df = pd.DataFrame(index=gene_names)
    adata = anndata.AnnData(X=results["pred"], obs=obs_df, var=var_df)
    # pred_de = torch.stack(pred_de)
    # truth_de = torch.stack(truth_de)
    # results["pred_de"] = pred_de.detach().cpu().numpy().astype(np.float32)
    # results["truth_de"] = truth_de.detach().cpu().numpy().astype(np.float32)

    return adata

def find_latest_checkpoint_dir(data_name: str, split_type: str, split_num: int) -> Path:
    """Finds the latest checkpoint directory for a given data_name and split_num."""
    # Pattern to match the saved model directories
    pattern = f"./save/{data_name}_{split_type}/dev_perturb_{data_name}-{split_num}-*"
    matching_dirs = glob.glob(pattern)

    if not matching_dirs:
        raise FileNotFoundError(f"No checkpoint directory found for {data_name} {split_type} split {split_num} matching pattern '{pattern}'")

    # Sort by modification time to get the latest
    matching_dirs.sort(key=os.path.getmtime, reverse=True)
    return Path(matching_dirs[0])

def main():
    parser = argparse.ArgumentParser(description="Predict expression profiles using trained scGPT models.")

    parser.add_argument(
        '--data_name',
        type=str,
        required=True,
        help='The name of the data (e.g., "k562", "rpe1") for which models were trained.'
    )
    parser.add_argument(
        '--split_type',
        type=str,
        required=True,
        help='The type of split (e.g. 10th, 90th)'
    )
    parser.add_argument(
        '--runs',
        type=str,
        required=True,
        help='Comma-separated list of run numbers (integers) to predict for (e.g., "0,1,2,4").'
    )
    parser.add_argument(
        '--eval_batch_size',
        type=int,
        default=64,
        help='Batch size for evaluation.'
    )
    parser.add_argument(
        '--pretrained_model_path',
        type=str, 
        default=os.path.join(PATH, "scripts/scGPT/save/scGPT_human"),
        help='Path to the directory containing the pretrained scGPT model (args.json, best_model.pt, vocab.json).'
    )
    
    args = parser.parse_args()
    
    # Directory to save prediction results
    prediction_results_dir = Path(f"/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/replogle_2022/{args.data_name}_{args.split_type}")
    prediction_results_dir.mkdir(parents=True, exist_ok=True)
    # Setting up logger
    logger = scg.logger
    scg.utils.add_file_handler(logger, prediction_results_dir / "run.log")
    logger.info(f"Saving prediction results to {prediction_results_dir}")
    
    finished_runs = [int(r.strip()) for r in args.runs.split(',')]
    print(f"Predicting for data: {args.data_name}, {args.split_type}, runs: {finished_runs}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    # Settings from training script that are needed for model initialization
    pad_token = "<pad>"
    special_tokens = [pad_token, "<cls>", "<eoc>"]
    pad_value = 0
    pert_pad_id = 0
    
    batch_size = args.eval_batch_size

    # Load vocab from the pretrained scGPT model (common for all runs)
    vocab_file = Path(args.pretrained_model_path) / "vocab.json"
    vocab = GeneVocab.from_file(vocab_file)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    logger.info(f"Loaded vocabulary of size {len(vocab)}")

    for run_num in finished_runs:
        logger.info(f"\n--- Processing run {run_num} for {args.data_name} {args.split_type} ---")
        
        try:
            model_save_dir = find_latest_checkpoint_dir(args.data_name, args.split_type, run_num)
            model_file = model_save_dir / "best_model.pt"
            
            # Use the model config from the pretrained model path
            model_config_file = Path(args.pretrained_model_path) / "args.json"
            if not model_config_file.exists():
                raise FileNotFoundError(f"Model config file not found: {model_config_file}")

            with open(model_config_file, "r") as f:
                model_configs = json.load(f)
            
            # Extract necessary model parameters from the loaded config
            embsize = model_configs.get("embsize", 512)
            nhead = model_configs.get("nheads", 8)
            d_hid = model_configs.get("d_hid", 512)
            nlayers = model_configs.get("nlayers", 12)
            n_layers_cls = model_configs.get("n_layers_cls", 3)
            dropout = model_configs.get("dropout", 0.0) # Assume default if not in config

            logger.info(f"Loading model from {model_file} with config from {model_config_file}")

            # Initialize model architecture
            ntokens = len(vocab)
            model = TransformerGenerator(
                ntokens,
                embsize,
                nhead,
                d_hid,
                nlayers,
                nlayers_cls=n_layers_cls,
                n_cls=1,
                vocab=vocab,
                dropout=dropout,
                pad_token=pad_token,
                pad_value=pad_value,
                pert_pad_id=pert_pad_id,
                use_fast_transformer=True, # Assuming this was used during training
            )
            
            model.load_state_dict(torch.load(model_file, map_location=device))
            model.to(device)
            logger.info(f"Model for run {run_num} loaded successfully.")

            # --- Data Loading for Evaluation ---
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
            pert_data_source_dataloader = pert_data_source.get_dataloader(batch_size=batch_size, test_batch_size=batch_size)
            pert_data_target.prepare_split(split = "no_split", seed = 42) 
            pert_data_target_dataloader=pert_data_target.get_dataloader(batch_size=batch_size, test_batch_size=batch_size)

            # Concatenate the datasets (target control + 1/10 or 9/10 of target unique perturbations)
            splits_df = pd.read_csv(os.path.join(PATH, "data/sigs/perturb-seq/pb_splits.csv")) # splits df is all in rpe1_shared and k562_shared
            split_col = f"split_{run_num}"

            # Get unique perturbations from the target dataset
            unique_perturbations = np.unique([data["pert_idx"][0] for data in pert_data_target_dataloader["test_loader"].dataset])

            # Extract perturbation names and idx from target
            pert_idx_to_name = {data["pert_idx"][0]: data["pert"].replace("+ctrl","") for data in pert_data_target_dataloader["test_loader"].dataset}
            pert_name_to_idx = {v: k for k, v in pert_idx_to_name.items()}

            # Get the split column name
            split_col = f"split_{run_num}"

            # Get perturbations in the eval split (TRUE values = held out 1/10th)
            if args.split_type == "90th":
                train_pb_names = splits_df[splits_df[split_col] == False]["pb"].values
                eval_pb_names = splits_df[splits_df[split_col] == True]["pb"].values
            elif args.split_type == "10th":
                train_pb_names = splits_df[splits_df[split_col] == True]["pb"].values
                eval_pb_names = splits_df[splits_df[split_col] == False]["pb"].values

            print(f"Using split {run_num}: {len(eval_pb_names)} eval perturbations, {len(train_pb_names)} train perturbations")

            # Convert names to indices (filter to only those present in the dataset)
            eval_perturbations = [pert_name_to_idx[name] for name in eval_pb_names if name in pert_name_to_idx]
            train_perturbations = [pert_name_to_idx[name] for name in train_pb_names if name in pert_name_to_idx]

            print(f"Found {len(eval_perturbations)} eval and {len(train_perturbations)} train perturbations in dataset")

            # Split the target dataset based on predefined splits
            train_data_split = [data for data in pert_data_target_dataloader["test_loader"].dataset 
                                if data["pert_idx"][0] in train_perturbations]
            eval_data_split = [data for data in pert_data_target_dataloader["test_loader"].dataset 
                               if data["pert_idx"][0] in eval_perturbations]

            combined_dataset = ConcatDataset([pert_data_source_dataloader["test_loader"].dataset, train_data_split])

            # Create a single DataLoader from the combined dataset
            combined_train_dataloader = DataLoader(combined_dataset, batch_size=batch_size, shuffle=True)
            valid_loader = DataLoader(eval_data_split, batch_size=batch_size, shuffle=True)
            
            logger.info(f"Loaded validation data for run {run_num} (length: {len(valid_loader.dataset)}).")
            
            # Derive gene_ids for eval_perturb (from source data and vocab, as in training)
            genes = pert_data_source.adata.var["gene_name"].tolist()
            gene_ids = np.array(
                [vocab[gene] if gene in vocab else vocab["<pad>"] for gene in genes], dtype=int
            )
            
            eval_configs = {
                "include_zero_gene": "all", # Assumed from training config
                "gene_ids": gene_ids,
                "gene_names": genes,
                "pert_list": pert_data_source.pert_names, # List of perts defined by set_pert_genes, pert_idx indexes into this
                "gene_map": pert_data_source.node_map # Dictionary mapping gene names to gene_var index in adata
            }

            # Run evaluation
            logger.info(f"Running eval_perturb for run {run_num}...")
            # results = eval_perturb(valid_loader, model, device, eval_configs, logger)
            pred_adata = eval_perturb(valid_loader, model, device, eval_configs, logger)
            logger.info(f"Prediction complete for run {run_num}.")

            # # Save results
            # result_filename = prediction_results_dir / f"predictions_{args.data_name}_run{run_num}.pt"
            # torch.save(results, result_filename)
            # logger.info(f"Saved predictions for run {run_num} to {result_filename}")
            
            # Save AnnData object
            adata_output_filename = prediction_results_dir / f"predicted_adata_{args.data_name}_{args.split_type}_run{run_num}.h5ad"
            pred_adata.write(adata_output_filename)
            logger.info(f"Saved AnnData object for run {run_num} to {adata_output_filename}")

        except FileNotFoundError as e:
            logger.error(f"Error processing run {run_num}: {e}")
        except Exception as e:
            logger.error(f"An unexpected error occurred for run {run_num}: {e}", exc_info=True)

    logger.info("\n--- All predictions complete ---")

if __name__ == "__main__":
    main()