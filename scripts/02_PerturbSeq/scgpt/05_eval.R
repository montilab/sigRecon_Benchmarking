library(tidyverse)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(PATH, "data/scgpt/prediction_results/")

# Loading data
rpe1_true_sigs <- sigrecon::perturbseq.rpe1
k562_true_sigs <- sigrecon::perturbseq.k562
stopifnot(all.equal(names(rpe1_true_sigs), names(k562_true_sigs)))

#10th sigs
rpe1_k562_pred_sigs <- readRDS(file.path(DATA_PATH, "rpe1_10th_sigs.rds"))
k562_rpe1_pred_sigs <- readRDS(file.path(DATA_PATH, "k562_10th_sigs.rds"))

names(rpe1_k562_pred_sigs) <- str_remove_all(names(rpe1_k562_pred_sigs), "\\+ctrl")
names(k562_rpe1_pred_sigs) <- str_remove_all(names(k562_rpe1_pred_sigs), "\\+ctrl")
rpe1_k562_filter <- intersect(names(rpe1_k562_pred_sigs), names(rpe1_true_sigs))
k562_rpe1_filter <- intersect(names(k562_rpe1_pred_sigs), names(k562_true_sigs))

rpe1_k562_pred_sigs <- rpe1_k562_pred_sigs[rpe1_k562_filter]
k562_rpe1_pred_sigs <- k562_rpe1_pred_sigs[k562_rpe1_filter]

# Eval tables
rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(rpe1_k562_pred_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs, 
                                       source = "rpe1_k562")
k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_rpe1_pred_sigs, function(x) x$up), 
                                       true_sigs = rpe1_true_sigs, 
                                       source = "k562_rpe1")

# No Change
rpe1_k562_eval_table_n <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         true_sigs = k562_true_sigs, 
                                         source = "rpe1_k562")
k562_rpe1_eval_table_n <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         true_sigs = rpe1_true_sigs, 
                                         source = "k562_rpe1")

all_eval_table <- rbind(rpe1_k562_eval_table,
                        k562_rpe1_eval_table)
all_eval_no_change_table <- rbind(rpe1_k562_eval_table_n,
                                  k562_rpe1_eval_table_n)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df, file.path(DATA_PATH, "10th_eval_table.rds"))

# Unpaired 
combined_aggregated_df %>% 
  dplyr::group_by(source) %>%
  dplyr::summarise(NES_mean = mean(NES_TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))

# Plots
combined_aggregated_df %>% 
  filter(source=="k562_rpe1") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "K562 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="k562_rpe1") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "K562 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1_k562") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "RPE1 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1_k562") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "RPE1 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")


# Paired Tests
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


#90th 
rpe1_k562_pred_sigs <- readRDS(file.path(DATA_PATH, "rpe1_90th_sigs.rds"))
k562_rpe1_pred_sigs <- readRDS(file.path(DATA_PATH, "k562_90th_sigs.rds"))

names(rpe1_k562_pred_sigs) <- str_remove_all(names(rpe1_k562_pred_sigs), "\\+ctrl")
names(k562_rpe1_pred_sigs) <- str_remove_all(names(k562_rpe1_pred_sigs), "\\+ctrl")
rpe1_k562_filter <- intersect(names(rpe1_k562_pred_sigs), names(rpe1_true_sigs))
k562_rpe1_filter <- intersect(names(k562_rpe1_pred_sigs), names(k562_true_sigs))

rpe1_k562_pred_sigs <- rpe1_k562_pred_sigs[rpe1_k562_filter]
k562_rpe1_pred_sigs <- k562_rpe1_pred_sigs[k562_rpe1_filter]

# Eval tables
rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(rpe1_k562_pred_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs, 
                                       source = "rpe1_k562")
k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_rpe1_pred_sigs, function(x) x$up), 
                                       true_sigs = rpe1_true_sigs, 
                                       source = "k562_rpe1")

# No Change
rpe1_k562_eval_table_n <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         true_sigs = k562_true_sigs, 
                                         source = "rpe1_k562")
k562_rpe1_eval_table_n <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         true_sigs = rpe1_true_sigs, 
                                         source = "k562_rpe1")

all_eval_table <- rbind(rpe1_k562_eval_table,
                        k562_rpe1_eval_table)
all_eval_no_change_table <- rbind(rpe1_k562_eval_table_n,
                                  k562_rpe1_eval_table_n)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df, file.path(DATA_PATH, "90th_eval_table.rds"))

# Unpaired 
combined_aggregated_df %>% 
  dplyr::group_by(source) %>%
  dplyr::summarise(NES_mean = mean(NES_TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))

# Plots
combined_aggregated_df %>% 
  filter(source=="k562_rpe1") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "K562 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="k562_rpe1") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "K562 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1_k562") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "RPE1 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1_k562") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "RPE1 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")


# Paired Tests
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))

