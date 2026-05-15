
load("C:\\Users\\abc\\Downloads\\yoga_data.RData")

library(ade4)
library(dplyr)
library(vegan)
library(pcaPP)
library(ggplot2)
library(ggrepel)

################################################################
## ALPHA DIVERSITY 
################################################################

shannon_index <- diversity(species_prof_tss, index = "shannon")

diversity_df <- as.data.frame(matrix(0, nrow=78, ncol=3))
rownames(diversity_df) <- rownames(species_prof_tss)
colnames(diversity_df) <- c("shannon","category","time_point")
diversity_df$category <- c(metadata_yoga_final$category,metadata_yoga_final$category)
diversity_df$time_point <- ifelse(1:nrow(diversity_df) <= 39, "baseline", "follow_up")
diversity_df$shannon <- shannon_index
diversity_df$adherence_groups <- c(metadata_yoga_final$adherence_groups,metadata_yoga_final$adherence_groups)

diversity_df$adherence_groups_new <- ifelse(diversity_df$category == 'walking','walking',ifelse(diversity_df$category == 'yoga' & diversity_df$adherence_groups == 1,'yoga1',ifelse(diversity_df$category == 'yoga' & diversity_df$adherence_groups == 2,"yoga2","yoga3")))

table(diversity_df$adherence_groups_new[1:39])

boxplot(diversity_df[(diversity_df$adherence_groups_new == 'walking' | diversity_df$adherence_groups_new == 'yoga1') & diversity_df$time_point == 'baseline','shannon'],
        diversity_df[(diversity_df$adherence_groups_new == 'walking' | diversity_df$adherence_groups_new == 'yoga1') & diversity_df$time_point == 'follow_up','shannon'],
        names = c("Baseline","Follow-up"), 
        col = c("skyblue1","gold1"), outline = F, ylab = "Shannon index")

wilcox.test(diversity_df[(diversity_df$adherence_groups_new == 'walking' | diversity_df$adherence_groups_new == 'yoga1') & diversity_df$time_point == 'baseline','shannon'],
            diversity_df[(diversity_df$adherence_groups_new == 'walking' | diversity_df$adherence_groups_new == 'yoga1') & diversity_df$time_point == 'follow_up','shannon'], paired = T)

boxplot(diversity_df[(diversity_df$adherence_groups_new == 'yoga2' | diversity_df$adherence_groups_new == 'yoga3') & diversity_df$time_point == 'baseline','shannon'],
        diversity_df[(diversity_df$adherence_groups_new == 'yoga2' | diversity_df$adherence_groups_new == 'yoga3') & diversity_df$time_point == 'follow_up','shannon'],
        names = c("Baseline","Follow-up"), 
        col = c("skyblue1","gold1"), outline = F, ylab = "Shannon index")

wilcox.test(diversity_df[(diversity_df$adherence_groups_new == 'yoga2' | diversity_df$adherence_groups_new == 'yoga3') & diversity_df$time_point == 'baseline','shannon'],
            diversity_df[(diversity_df$adherence_groups_new == 'yoga2' | diversity_df$adherence_groups_new == 'yoga3') & diversity_df$time_point == 'follow_up','shannon'])

#################################################################
# Association with Simple Clinical Colitis Activity Index (SCCAI)
#################################################################

# Select only samples ending with "_A"
yoga_samples <- grepl("A_", rownames(species_prof_tss))

species_prof_yoga <- species_prof_tss[yoga_samples, ]
species_prof_yoga_bs <- species_prof_yoga[grepl("_B$", rownames(species_prof_yoga)), ]

# Bray-Curtis distance
bray_dist_yoga_bs <- vegdist(species_prof_yoga_bs, method = "bray")

yoga_bs_pcoa <- dudi.pco(bray_dist_yoga_bs, scannf = F, nf = 2)

yoga_bs_pcoa_df <- as.data.frame(yoga_bs_pcoa$li)
colnames(yoga_bs_pcoa_df) <- c("PCoA1", "PCoA2")

yoga_bs_pcoa_df$SCCAI <- metadata_yoga_final[metadata_yoga_final$category == 'yoga',' SCCAI']

permanova_result <- adonis2(
  yoga_bs_pcoa_df[, c("PCoA1", "PCoA2")] ~ SCCAI,
  data = yoga_bs_pcoa_df,
  method = "euclidean",
  permutations = 999)

permanova_result


############# After Week 24 #############

species_prof_yoga_w24 <- species_prof_yoga[grepl("_F$", rownames(species_prof_yoga)), ]

# Bray-Curtis distance
bray_dist_yoga_w24 <- vegdist(species_prof_yoga_w24, method = "bray")

yoga_w24_pcoa <- dudi.pco(bray_dist_yoga_w24, scannf = F, nf = 2)

yoga_w24_pcoa_df <- as.data.frame(yoga_w24_pcoa$li)
colnames(yoga_w24_pcoa_df) <- c("PCoA1", "PCoA2")

yoga_w24_pcoa_df$SCCAI <- metadata_yoga_final[metadata_yoga_final$category == 'yoga','24 WEEKKS SCCAI']

permanova_result <- adonis2(
  yoga_w24_pcoa_df[, c("PCoA1", "PCoA2")] ~ SCCAI,
  data = yoga_w24_pcoa_df,
  method = "euclidean",
  permutations = 999)

permanova_result

##At baseline, patients likely had stronger disease-associated dysbiosis
##so microbiome variation tracked SCCAI severity


##################################################################
## volcano-plot of taxa that show significant changes in walking + yoga1 between week24 and baseline
##################################################################

walking_yoga1_SpeciesProfile_B <- species_prof_tss[c(rownames(diversity_df[diversity_df$adherence_groups_new == 'walking' & diversity_df$time_point == 'baseline', ]),rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga1' & diversity_df$time_point == 'baseline',])),]
walking_yoga1_SpeciesProfile_F <- species_prof_tss[c(rownames(diversity_df[diversity_df$adherence_groups_new == 'walking' & diversity_df$time_point == 'follow_up', ]),rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga1' & diversity_df$time_point == 'follow_up',])),]

wilcox_batch_paired_mean = function(x,y)
{
  p_array <- NULL;
  type_array <- NULL;
  mean1_array <- NULL;
  mean2_array <- NULL;
  x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
  y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
  z <- intersect(rownames(x),rownames(y));
  for(i in 1:length(z))
  {
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]),paired=TRUE)$p.value;
    type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
    mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
    i <- i + 1;
  }
  out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array, method = 'BH'),mean1_array,mean2_array));
  rownames(out) <- z;
  colnames(out)[3] <- "p_adj"
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  return(t(out));
}

wilcox_paired_mean_walking_yoga1_F_B <- as.data.frame(wilcox_batch_paired_mean(t(walking_yoga1_SpeciesProfile_F),t(walking_yoga1_SpeciesProfile_B)))

wilcox_paired_mean_walking_yoga1_F_B$log2_fold_change <- log2(wilcox_paired_mean_walking_yoga1_F_B$mean1_array / wilcox_paired_mean_walking_yoga1_F_B$mean2_array)
wilcox_paired_mean_walking_yoga1_F_B$transformed_p <- -log10(wilcox_paired_mean_walking_yoga1_F_B$p_array)
wilcox_paired_mean_walking_yoga1_F_B$color_code <- ifelse(wilcox_paired_mean_walking_yoga1_F_B$log2_fold_change > 1, "color1",ifelse(wilcox_paired_mean_walking_yoga1_F_B$log2_fold_change < -1, "color2","color3"))
wilcox_paired_mean_walking_yoga1_F_B_trimmed <- wilcox_paired_mean_walking_yoga1_F_B[wilcox_paired_mean_walking_yoga1_F_B$transformed_p >1,]

ggplot(wilcox_paired_mean_walking_yoga1_F_B_trimmed, aes(x = log2_fold_change, y = transformed_p)) +
  geom_point(aes(color = color_code), size = 3) +
  geom_text_repel(
    data = subset(wilcox_paired_mean_walking_yoga1_F_B_trimmed, color_code %in% c("color1", "color2")),
    aes(label = rownames(subset(wilcox_paired_mean_walking_yoga1_F_B_trimmed, color_code %in% c("color1", "color2")))),
    size = 3.5,
    max.overlaps = 45
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  scale_color_manual(values = c("color1" = "royalblue3", "color2" = "deeppink1", "color3" = "lightgray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (Follow-up/Baseline)", y = "-Log10 (p-value)") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14))


##################################################################
## volcano-plot of taxa that show significant changes in yoga2 + yoga3 between week24 and baseline
##################################################################

yoga2_yoga3_SpeciesProfile_B <- species_prof_tss[c(rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga2' & diversity_df$time_point == 'baseline', ]),rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga3' & diversity_df$time_point == 'baseline',])),]
yoga2_yoga3_SpeciesProfile_F <- species_prof_tss[c(rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga2' & diversity_df$time_point == 'follow_up', ]),rownames(diversity_df[diversity_df$adherence_groups_new == 'yoga3' & diversity_df$time_point == 'follow_up',])),]

wilcox_paired_mean_yoga2_yoga3_F_B <- as.data.frame(wilcox_batch_paired_mean(t(yoga2_yoga3_SpeciesProfile_F),t(yoga2_yoga3_SpeciesProfile_B)))

wilcox_paired_mean_yoga2_yoga3_F_B$log2_fold_change <- log2(wilcox_paired_mean_yoga2_yoga3_F_B$mean1_array / wilcox_paired_mean_yoga2_yoga3_F_B$mean2_array)
wilcox_paired_mean_yoga2_yoga3_F_B$transformed_p <- -log10(wilcox_paired_mean_yoga2_yoga3_F_B$p_array)
wilcox_paired_mean_yoga2_yoga3_F_B$color_code <- ifelse(wilcox_paired_mean_yoga2_yoga3_F_B$log2_fold_change > 1, "color1",ifelse(wilcox_paired_mean_yoga2_yoga3_F_B$log2_fold_change < -1, "color2","color3"))
wilcox_paired_mean_yoga2_yoga3_F_B_trimmed <- wilcox_paired_mean_yoga2_yoga3_F_B[wilcox_paired_mean_yoga2_yoga3_F_B$transformed_p >1,]

ggplot(wilcox_paired_mean_yoga2_yoga3_F_B_trimmed, aes(x = log2_fold_change, y = transformed_p)) +
  geom_point(aes(color = color_code), size = 3) +
  geom_text_repel(
    data = subset(wilcox_paired_mean_yoga2_yoga3_F_B_trimmed, color_code %in% c("color1", "color2")),
    aes(label = rownames(subset(wilcox_paired_mean_yoga2_yoga3_F_B_trimmed, color_code %in% c("color1", "color2")))),
    size = 3.5,
    max.overlaps = 45
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  scale_color_manual(values = c("color1" = "royalblue3", "color2" = "deeppink1", "color3" = "lightgray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (Follow-up/Baseline)", y = "-Log10 (p-value)") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14))

################################################################
## KENDALL Correlation 
################################################################

Yoga_metadata_F <- subset(metadata_yoga_final[metadata_yoga_final$category == 'yoga',], select = c(adherence))
Yoga_metadata_F$sample_id <- rownames(species_prof_yoga_w24)
rownames(Yoga_metadata_F) <- Yoga_metadata_F$sample_id

####### spearman correlation
Yoga_species_profile_F <- species_prof_tss[rownames(Yoga_metadata_F),]
Yoga_species_profile_F <- Yoga_species_profile_F[,which(colSums(Yoga_species_profile_F)!=0)]
Yoga_species_profile_F <- Yoga_species_profile_F[which(rowSums(Yoga_species_profile_F)!=0),]

Yoga_species_profile_F$sample_id <- rownames(Yoga_species_profile_F)
Yoga_species_profile_F <- left_join(Yoga_species_profile_F,Yoga_metadata_F,by = "sample_id")
rownames(Yoga_species_profile_F) <- Yoga_species_profile_F$sample_id
Yoga_species_profile_F$sample_id <- NULL


## create an empty dataframe to store the correlation results for spearman
corr_F_sp_Adherence <- data.frame(Species = character(), Corr_estimate = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

for (species in colnames(Yoga_species_profile_F)[colnames(Yoga_species_profile_F) != "adherence"]) {
  cor_test <- cor.test(Yoga_species_profile_F[[species]], Yoga_species_profile_F$adherence, method = "spearman", exact = NULL)
  corr_F_sp_Adherence <- rbind(corr_F_sp_Adherence, data.frame(Species = species, Corr_estimate = cor_test$estimate, PValue = cor_test$p.value))
}


# Add a column for -log10(p-value)
corr_F_sp_Adherence$logPValue <- -log10(corr_F_sp_Adherence$PValue)
corr_F_sp_Adherence$Significant <- ifelse(corr_F_sp_Adherence$PValue <= 0.05, "significant_pval", "non-significant_pval")



volcano_plot <- ggplot(corr_F_sp_Adherence, aes(x = Corr_estimate, y = logPValue)) +
  geom_point(aes(color = Significant)) +
  geom_text_repel(data = subset(corr_F_sp_Adherence, Significant == "significant_pval"),
                  aes(label = Species), size = 5, max.overlaps = 45) +  # Adjust direction based on Corr_estimate
  scale_color_manual(values = c("non-significant_pval" = "black", "significant_pval" = "red")) +
  theme_minimal() +
  labs(title = "Correlation of Species at F with Adherence", x = "Spearman Corr_estimate", y = "-log10(p-value)") +
  theme(panel.border = element_rect(color = "black", fill = NA,size = 2))

volcano_plot


###################################################################
## Correlating microbiome Loss with Inflammation Markers using 
## follow-up - baseline values 
##################################################################

commonTaxa_2groups_paired_wilcox <- intersect(rownames(wilcox_paired_mean_walking_yoga1_F_B), rownames(wilcox_paired_mean_yoga2_yoga3_F_B))

df_inflammatory_markers_diff <- data.frame(category = metadata_yoga_final$category,
                                           ALB = metadata_yoga_final$ALB2 - metadata_yoga_final$ALB,
                                           CRP = metadata_yoga_final$CRP2 - metadata_yoga_final$CRP,
                                           ESR = metadata_yoga_final$ESR2 - metadata_yoga_final$ESR, 
                                           FCP = as.numeric(metadata_yoga_final$`FCP F`) - as.numeric(metadata_yoga_final$`FCP B`), 
                                           PSS = metadata_yoga_final$`PSS-2` - metadata_yoga_final$PSS, 
                                           HAM = metadata_yoga_final$`HAM-2` - metadata_yoga_final$`HAM-1`, 
                                           adherence_groups = metadata_yoga_final$adherence_groups,
                                           row.names = rownames(metadata_yoga_final))

species_prof_tss_selected_diff <- species_prof_tss[40:78,commonTaxa_2groups_paired_wilcox] - species_prof_tss[1:39,commonTaxa_2groups_paired_wilcox]


species_mat <- as.data.frame(species_prof_tss_selected_diff)  
marker_mat <- df_inflammatory_markers_diff[, 2:7] 
marker_mat[] <- lapply(marker_mat, function(x) as.numeric(as.character(x)))

corr_mat <- matrix(NA, nrow = ncol(species_mat), ncol = ncol(marker_mat))
pval_mat <- matrix(NA, nrow = ncol(species_mat), ncol = ncol(marker_mat))
dir_mat  <- matrix(0, nrow = ncol(species_mat), ncol = ncol(marker_mat))

rownames(corr_mat) <- rownames(pval_mat) <- rownames(dir_mat) <- colnames(species_mat)
colnames(corr_mat) <- colnames(pval_mat) <- colnames(dir_mat) <- colnames(marker_mat)

for (i in 1:ncol(species_mat)) {
  for (j in 1:ncol(marker_mat)) {
    sp_values <- species_mat[[i]]
    mk_values <- marker_mat[[j]]
    
    valid_idx <- which(!is.na(sp_values) & !is.na(mk_values))
    
    if (length(valid_idx) > 2) {  
      ct <- cor.test(sp_values[valid_idx], mk_values[valid_idx], method = "spearman")
      corr_mat[i, j] <- ct$estimate
      pval_mat[i, j] <- ct$p.value
      
      if (!is.na(ct$p.value)) {
        if (ct$p.value <= 0.05) {
          dir_mat[i, j] <- ifelse(ct$estimate > 0, 2, -2)
        } else if (ct$p.value <= 0.1) {
          dir_mat[i, j] <- ifelse(ct$estimate > 0, 1, -1)
        }
      }
    }
  }
}

corr_df <- as.data.frame(corr_mat)
pval_df <- as.data.frame(pval_mat)
dir_df  <- as.data.frame(dir_mat)

inflammatory_markers_corr_results_diff <- list(correlation = corr_df, pvalue = pval_df, direction = dir_df)

dir_df_filtered_diff <- inflammatory_markers_corr_results_diff[["direction"]][rowSums(inflammatory_markers_corr_results_diff[["direction"]] != 0, na.rm = TRUE) >= 2,]

##### Making the heatmap ####

dir_df_filtered_diff_t <- t(dir_df_filtered_diff)

colnames(dir_df_filtered_diff_t) <- gsub("_", " ", colnames(dir_df_filtered_diff_t))

dir_df_filtered_t_num_diff <- apply(dir_df_filtered_diff_t, 2, function(x) as.numeric(as.character(x)))
rownames(dir_df_filtered_t_num_diff) <- rownames(dir_df_filtered_diff_t)

custom_colors <- c("#377eb8", "#a6cee3", "white", "#fb9a99", "#e31a1c")  # 5 colors

breaks <- seq(-2.5, 2.5, length.out = 6)

pheatmap::pheatmap(
  dir_df_filtered_t_num_diff,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,       
  treeheight_col = 0,       
  color = custom_colors,
  breaks = breaks,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellheight = 10,
  cellwidth = 10,
  main = "Associating taxa with Inflammatory markers")

# F.prausnitzii (FCP & ESR neg), D.longicatena (ESR neg) and Lachnospira eligens(PSS neg)
# M. faecis (PSS neg), F. saccharivorans (ESR neg)








