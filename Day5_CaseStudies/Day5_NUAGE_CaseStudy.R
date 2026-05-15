load("NuAGE_data.RData")

###  load packages  ###
library(vegan)
library(ade4)
library(pcaPP)
library(randomForest)
library(ggplot2)
library(dplyr)

# ─── STEP 1.1: Get to know your data ─────────────────────────────────────────
# Always start by understanding the dimensions and structure of your tables.

# How many samples and features do we have?
cat("=== DATASET DIMENSIONS ===\n")
cat("Metadata:         ", nrow(NuAGE_Metadata), "samples x", ncol(NuAGE_Metadata), "variables\n")
cat("Species profile:  ", nrow(NuAGE_SpProfile), "samples x", ncol(NuAGE_SpProfile), "species\n")
cat("Food intake:      ", nrow(NuAGE_FoodIntakeProfile), "samples x", ncol(NuAGE_FoodIntakeProfile), "nutrients\n")

# Peek at the metadata
head(NuAGE_Metadata)

# What groups do we have?
table(NuAGE_Metadata$study_condition)

# What countries are represented?
table(NuAGE_Metadata$country)

# Check age categories
table(NuAGE_Metadata$age_category)

# The NU-AGE study is LONGITUDINAL: each person was measured BEFORE and AFTER
# the 1-year dietary intervention. We need to track this in our data.
table(sub(".*_", "", rownames(NuAGE_Metadata))) #T0: 610; T1: 610
#nothing to drop

# ─── STEP 1.2: Split the data based on BEFORE and AFTER ─────────────────────────────────────────
NuAGE_SpProfile_Before <- NuAGE_SpProfile[grep("_T0$", rownames(NuAGE_SpProfile)), ]
NuAGE_SpProfile_After <- NuAGE_SpProfile[grep("_T1$", rownames(NuAGE_SpProfile)), ]

rownames(NuAGE_SpProfile_Before) <- gsub("_T0","", rownames(NuAGE_SpProfile_Before))
rownames(NuAGE_SpProfile_After) <- gsub("_T1","", rownames(NuAGE_SpProfile_After))

NuAGE_Metadata_Before <- NuAGE_Metadata[grep("_T0$", rownames(NuAGE_Metadata)), ]
NuAGE_Metadata_After <- NuAGE_Metadata[grep("_T1$", rownames(NuAGE_Metadata)), ]

rownames(NuAGE_Metadata_Before) <- gsub("_T0","", rownames(NuAGE_Metadata_Before))
rownames(NuAGE_Metadata_After) <- gsub("_T1","", rownames(NuAGE_Metadata_After))

NuAGE_FoodIntakeProfile_Before <- NuAGE_FoodIntakeProfile[grep("_T0$", rownames(NuAGE_FoodIntakeProfile)), ]
NuAGE_FoodIntakeProfile_After <- NuAGE_FoodIntakeProfile[grep("_T1$", rownames(NuAGE_FoodIntakeProfile)), ]

rownames(NuAGE_FoodIntakeProfile_Before) <- gsub("_T0","", rownames(NuAGE_FoodIntakeProfile_Before))
rownames(NuAGE_FoodIntakeProfile_After) <- gsub("_T1","", rownames(NuAGE_FoodIntakeProfile_After))

# ─── STEP 2: ALPHA DIVERSITY  ─────────────────────────────────────────
### Computing Alpha diversity

## Baseline 16s species profile
shannon_index_16s <- diversity(NuAGE_SpProfile_Before, index = "shannon")

shannon_df_16s <- data.frame(
  Sample = rownames(NuAGE_SpProfile_Before),
  Shannon_Index = shannon_index_16s,
  country = NuAGE_Metadata_Before$country)

# Baseline food intake 
shannon_index_food <- diversity(NuAGE_FoodIntakeProfile_Before, index = "shannon")

shannon_df_food <- data.frame(
  Sample = rownames(NuAGE_FoodIntakeProfile_Before),
  Shannon_Index = shannon_index_food,
  country = NuAGE_Metadata_After$country)

boxplot(shannon_df_16s$Shannon_Index,shannon_df_food$Shannon_Index, names = c("16s","food_intake"), col = c("royalblue","gold"), outline = F)
wilcox.test(shannon_df_16s$Shannon_Index,shannon_df_food$Shannon_Index, paired = T)

# ─── STEP 3: BETA DIVERSITY  ─────────────────────────────────────────
### Between Baseline 16s species profile and Baseline food intake 

### 16S SPECIES PROFILE ###

mat_species <- as.matrix(NuAGE_SpProfile_Before)

dist_species <- as.dist(1 - cor.fk(t(mat_species)) / 2)

pcoa_species <- dudi.pco(dist_species, scannf = FALSE, nf = 20)

pcoa_species_df <- as.data.frame(pcoa_species$li[,1:2])

colnames(pcoa_species_df) <- c("PCoA1", "PCoA2")

pcoa_species_df$country <- NuAGE_Metadata_Before$country

grp <- as.factor(pcoa_species_df$country)

s.class(
  pcoa_species_df,
  fac = grp,
  col = c("red","blue","green","orange","purple"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE,
  cpoint  = 2
)

adonis_species <- adonis2(dist_species ~ country,
                          data = pcoa_species_df)

adonis_species

### FOOD INTAKE PROFILE ###

mat_food <- as.matrix(NuAGE_FoodIntakeProfile_Before)

dist_food <- as.dist(1 - cor.fk(t(mat_food)) / 2)

pcoa_food <- dudi.pco(dist_food, scannf = FALSE, nf = 20)

pcoa_food_df <- as.data.frame(pcoa_food$li[,1:2])

colnames(pcoa_food_df) <- c("PCoA1", "PCoA2")

pcoa_food_df$country <- NuAGE_Metadata_Before$country

grp <- as.factor(pcoa_food_df$country)

s.class(
  pcoa_food_df,
  fac = grp,
  col = c("red","blue","green","orange","purple"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE,
  cpoint  = 2
)

adonis_food <- adonis2(dist_food ~ country,
                       data = pcoa_food_df)

adonis_food

### PROCRUSTES ###

procuste.randtest(pcoa_species$li[,1:2], pcoa_food$li[,1:2])

# Inference: The results indicate that there are country-specific patterns in 
# dietary habits which are also reflected in the microbiome profiles.

# ─── STEP 4: Random Forest Model  ─────────────────────────────────────────
## Identification of diet responsive taxa by machine learning

### RANDOM FOREST ###
x <- NuAGE_SpProfile
y <- NuAGE_Metadata$food_scores

RF_result <- randomForest(x, y)

# ACTUAL VS PREDICTED
predicted_scores <- RF_result$predicted

cor.test(predicted_scores, y, method = "spearman")

plot(y,
     predicted_scores,
     pch = 16,
     xlab = "Actual Diet Score",
     ylab = "Predicted Diet Score")

### FEATURE IMPORTANCE ###
rank_scale=function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

FI_scaled <- data.frame(
  Taxa = rownames(RF_result$importance),
  Overall = rank_scale(RF_result$importance[,1])
)

### SORT FEATURES ###
FI_scaled <- FI_scaled[order(FI_scaled$Overall, decreasing = TRUE), ]
head(FI_scaled)

### CUMULATIVE IMPORTANCE ###
FI_scaled$Percent <- FI_scaled$Overall / sum(FI_scaled$Overall)
FI_scaled$Cumulative <- cumsum(FI_scaled$Percent)

Top_species <- FI_scaled[FI_scaled$Cumulative <= 0.90, ]
head(Top_species)

### CORRELATION ###
cor_results <- apply(NuAGE_SpProfile[, Top_species$Taxa],2,
  function(x){
    test <- cor.test(x, y, method = "spearman")
    c(Rho = test$estimate,
      Pvalue = test$p.value)
  }
)

cor_results <- as.data.frame(t(cor_results))

Top_species$Rho <- cor_results$Rho
Top_species$Pvalue <- cor_results$Pvalue

Top_species_sig <- Top_species[Top_species$Pvalue <= 0.05, ]

### DIRECTION ###
Top_species_sig$Diet_group <- ifelse(Top_species_sig$Rho > 0,"Positive", "Negative")

### BARPLOT ###
# CREATE GROUPS
positive_df <- Top_species_sig %>% filter(Rho > 0) %>% arrange(Overall)

negative_df <- Top_species_sig %>% filter(Rho < 0) %>% arrange(desc(Overall))

# ADD LABEL
positive_df$Group <- "Diet Positive"
negative_df$Group <- "Diet Negative"

#showing just the top 25 species 
positive_df <- positive_df %>% slice_max(Overall, n = 25)
negative_df <- negative_df %>% slice_max(Overall, n = 25)
plot_df <- rbind(positive_df, negative_df)

# PLOT
ggplot(plot_df,
       aes(
         x = reorder(Taxa, Overall),
         y = Overall,
         fill = Group
       )) +
  
  geom_bar(stat = "identity") +
  
  coord_flip() +
  
  facet_wrap(~Group, scales = "free_y") +
  
  scale_fill_manual(values = c(
    "Diet Positive" = "steelblue",
    "Diet Negative" = "darkorange"
  )) +
  
  labs(
    title = "Top Important Species",
    x = "",
    y = "Overall Importance (Scaled)"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(size = 7),
    
    legend.position = "none"
  )


#correlation with AdherenceScores

### SELECT SPECIES ###
species_name <- "Faecalibacterium_prausnitzii" #positive
#species_name <- "Eubacterium_xylanophilum" #positive
#species_name <- "Ruminococcus_torques" #negative
#species_name <- "Collinsella_aerofaciens" #negative

### CREATE TEMP DATAFRAME ###
plot_df <- data.frame(
  
  Adherence =
    NuAGE_Metadata$AdherenceScores,
  
  Species =
    NuAGE_SpProfile[, species_name]
)

### SPEARMAN CORRELATION ###
test <- cor.test(
  plot_df$Adherence,
  plot_df$Species,
  method = "spearman"
)

### EXTRACT VALUES ###
rho <- round(test$estimate, 3)

pval <- signif(test$p.value, 3)

### SCATTERPLOT ###
ggplot(
  plot_df,
  aes(
    x = Adherence,
    y = Species
  )
) +
  
  geom_point(
    size = 2,
    alpha = 0.7,
    color = "steelblue"
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "darkred"
  ) +
  
  labs(
    title = species_name,
    
    subtitle = paste0(
      "Spearman rho = ",
      rho,
      " | p = ",
      pval
    ),
    
    x = "Mediterranean Diet Adherence Score",
    
    y = "Species Abundance"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )











