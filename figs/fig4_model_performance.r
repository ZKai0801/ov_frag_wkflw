library(tidyverse)
library(data.table)
library(readxl)
library(pROC)

getwd()
setwd("\\\\172.16.11.242/test_data/frag_OV/results/features")
wrkdir <- "\\\\172.16.11.242/test_data/frag_OV/manuscript/version3"

cancer_scores <- fread("models3/prediction_unenriched.csv")
groups <- as.data.table(read_excel(file.path(wrkdir, "/Supplementary_Tables.xlsx")))

setkey(groups, sampleID)
setkey(cancer_scores, sampleID)

df <- cancer_scores[groups]



# -----------------------------------
# Fig3A AUC

library(pROC)

cnv <- fread("./prediction_cnv.csv")
edm <- fread("./prediction_edm.csv")
bpm <- fread("./prediction_bpm.csv")
fsr <- fread("./prediction_fsr.csv")
fsd <- fread("./prediction_fsd.csv")


model_comb <- roc(df$group, df$proba_1)
model_cnv <- roc(cnv$group, cnv$proba_1)
model_edm <- roc(edm$group, edm$proba_1)
model_bpm <- roc(bpm$group, bpm$proba_1)
model_fsr <- roc(fsr$group, fsr$proba_1)
model_fsd <- roc(fsd$group, fsd$proba_1)


auc(df$group, df$proba_1)
auc(cnv$group, cnv$proba_1)
auc(edm$group, edm$proba_1)
auc(bpm$group, bpm$proba_1)
auc(fsr$group, fsr$proba_1)
auc(fsd$group, fsd$proba_1)


ci.auc(df$group, df$proba_1)
ci.auc(cnv$group, cnv$proba_1)
ci.auc(edm$group, edm$proba_1)
ci.auc(bpm$group, bpm$proba_1)
ci.auc(fsr$group, fsr$proba_1)
ci.auc(fsd$group, fsd$proba_1)




svg(file.path(wrkdir, "./test.svg"), 
    height = 5, width = 5)

plot(model_comb, col = "#4169E1", 
     xlab = "False Positive Rate", 
     ylab = "True Positive Rate",
     asp = NA,
     grid = TRUE,
     legacy.axes = TRUE)

plot(model_cnv, col = "#FF7F50", add = T)
plot(model_edm, col = "#20B2AA", add = T)
plot(model_bpm, col = "#9ACD32", add = T)
plot(model_fsr, col = "#FFA500", add = T)
plot(model_fsd, col = "#FA8072", add = T)


legend("bottomright",
       legend = c("Stacked:0.974 (95% CI: 0.967-0.980)", 
                  "CNV:0.901 (95% CI: 0.887-0.914)", 
                  "EDM:0.918 (95% CI: 0.905-0.931)", 
                  "BPM:0.919 (95% CI: 0.906-0.932)", 
                  "FSR:0.881 (95% CI: 0.867-0.896)", 
                  "FSD:0.923 (95% CI: 0.911-0.936)"),
       col = c("#4169E1", "#FF7F50", "#20B2AA", 
               "#9ACD32","#FFA500", "#FA8072"),
       lwd=4, cex =0.6, xpd = TRUE)

dev.off()


# -----------------------------------
# Fig3B

# figo I/II & III/IV

all_sens <- c()
for (irange in seq(1,10)){
    subdf <- df %>%
        filter(fold_index >= (irange-1)*10, 
               fold_index < irange*10) %>%
        arrange(desc(proba_1))
    
    sens_spec_df <- data.frame(threshold = numeric(), 
                               sensitivity = numeric(), 
                               specificity = numeric())
    
    for (i in seq(0, 1, 0.01)) {
        pred <- ifelse(subdf$proba_1 >= i, 1, 0)
        pred2 <- ifelse(subdf %>% filter(group ==1, (FIGO_stage == "I" | FIGO_stage == "II"))  %>% pull(proba_1) > i, 1, 0)

        if (sum(pred == 1) == 0 || sum(pred == 0) == 0) {
            cm <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(c("0", "1"), c("0", "1")))
        } else {
            cm <- table(pred, subdf$group)
        }
        
        specificity <- cm[1, 1] / sum(cm[1, ])
        sensitivity <- sum(pred2) / subdf %>% filter(group == 1, (FIGO_stage == "I" | FIGO_stage == "II"))  %>% nrow()
        sens_spec_df <- rbind(sens_spec_df, data.frame(threshold = i, sensitivity = sensitivity, specificity = specificity))
    }
    
    sens <- sens_spec_df %>%
        arrange(abs(specificity - 0.9)) %>%
        head(n = 1) %>%
        pull("sensitivity") 
    
    all_sens <- append(all_sens, sens)
}



figo_res %>%
    as_tibble() %>%
    ggplot(aes(x = figo, y = sens)) +
    geom_bar(aes(fill = figo), stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = sens-ci, ymax = sens+ci),
                  width = 0.1) +
    geom_text(aes(y = 0.06, label = round(sens,2))) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "FIGO stages", y = "Sensitivity") +
    scale_fill_manual(values = c("#fc8d59", "#e34a33"))


ggsave(file.path(wrkdir, "Fig3C_FIGO_comb.pdf"),
       width = 5, height = 5)



# -------------------------------------------
# Fig. 4C

df <- df %>%
    mutate(Diagnosis = case_when(grepl("Hgsc", Diagnosis) ~ "HGSC",
                                 grepl("Lgsc", Diagnosis) ~ "LGSC",
                                 grepl("Endometrioid", Diagnosis) ~ "Endometrioid",
                                 grepl("Clear Cell", Diagnosis) ~ "CCO",
                                 T ~ "Others"))


all_sens <- c()
for (irange in seq(1,10)){
    subdf <- df %>%
        filter(fold_index >= (irange-1)*10, 
               fold_index < irange*10) %>%
        arrange(desc(proba_1))
    
    sens_spec_df <- data.frame(threshold = numeric(), 
                               sensitivity = numeric(), 
                               specificity = numeric())
    
    for (i in seq(0, 1, 0.01)) {
        pred <- ifelse(subdf$proba_1 >= i, 1, 0)
        pred2 <- ifelse(subdf %>% filter(group ==1, Diagnosis == "Others")  %>% pull(proba_1) > i, 1, 0)
        
        if (sum(pred == 1) == 0 || sum(pred == 0) == 0) {
            cm <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(c("0", "1"), c("0", "1")))
        } else {
            cm <- table(pred, subdf$group)
        }
        
        specificity <- cm[1, 1] / sum(cm[1, ])
        sensitivity <- sum(pred2) / subdf %>% filter(group ==1, Diagnosis == "Others")  %>% nrow()
        sens_spec_df <- rbind(sens_spec_df, data.frame(threshold = i, sensitivity = sensitivity, specificity = specificity))
    }
    
    sens <- sens_spec_df %>%
        arrange(abs(specificity - 0.9)) %>%
        head(n = 1) %>%
        pull("sensitivity")    
    all_sens <- append(all_sens, sens)
}

histo_subtypes %>%
    as_tibble() %>%
    mutate(histo = factor(histo, level = c("HGSC", "LGSC", "CCO", "Endometroid", "Others"))) %>%
    ggplot(aes(x = histo, y = sens)) +
    geom_bar(aes(fill = histo), stat = "identity") +
    geom_errorbar(aes(ymin = sens-ci, ymax = sens+ci),
                  width = 0.2) +
    geom_text(aes(y = 0.06, label = round(sens,2))) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", y = "Sensitivity") +
    ggsci::scale_fill_npg()

ggsave(file.path(wrkdir, "Fig3C_Histo.png"),
       width = 5, height = 5)


# -----------------------------------
# Fig4d mean cancer scores

p_overall <- df %>%
    mutate(groups = case_when(i.group == "OV" ~ "OV",
                              i.group == "Benigns" ~ "Benign",
                              i.group == "Healths" ~ "Healthy")) %>%
    mutate(groups = factor(groups, level = c("OV", "Benign", "Healthy"))) %>%
    group_by(sampleID, groups) %>%
    summarise(mean_prob = mean(proba_1)) %>%
    ggplot(aes(x = groups, y = mean_prob,
               color = groups, fill = groups)) + 
    geom_violin(alpha = 0.7) +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_color_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    labs(x = "", y = "Mean Cancer Scores") +
    theme(legend.position = "none")


ggsave(file.path(wrkdir, "Fig3D_Mean_cancerscores.png"),
       width = 5, height = 5)



# -----------------------------------
# Fig 4E.

tfx <- fread("./ichorCNA_tfx.tsv")

options(scipen=999)

p_tfx <- df %>%
    mutate(groups = case_when(i.group == "OV" ~ "OV",
                              i.group == "Benigns" ~ "Benign",
                              i.group == "Healths" ~ "Healthy")) %>%
    mutate(groups = factor(groups, level = c("OV", "Benign", "Healthy"))) %>%
    group_by(sampleID, groups) %>%
    summarise(mean_proba = mean(proba_1)) %>%
    left_join(tfx[, list(sampleID, V2, V3)]) %>%
    filter(groups == "OV") %>%
    mutate(group2 = if_else(V2 > 0.03, "TFX > 0.03", "TFX <= 0.03")) %>%
    ggplot(aes(x = group2, y = mean_proba, fill = group2, color = group2)) +
    geom_violin(alpha = 0.7) +
    geom_point() +
    theme_classic() +
    scale_fill_manual(values = c("#4682B4", "#DC143C")) +
    scale_color_manual(values = c("#4682B4", "#DC143C")) +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    ggpubr::stat_compare_means(comparisons = list(c("TFX > 0.03", "TFX <= 0.03")),
                               label = "p.format")  +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


p_figo <- df %>%
    mutate(groups = case_when(i.group == "OV" ~ "OV",
                              i.group == "Benigns" ~ "Benign",
                              i.group == "Healths" ~ "Healthy")) %>%
    mutate(groups = factor(groups, level = c("OV", "Benign", "Healthy"))) %>%
    group_by(sampleID, groups, FIGO_stage) %>%
    summarise(mean_proba = mean(proba_1)) %>%
    mutate(FIGO_stage = case_when(FIGO_stage == "I" ~ "I/II",
                                  FIGO_stage == "II" ~ "I/II",
                                  FIGO_stage == "III" ~ "III/IV",
                                  FIGO_stage == "IV" ~ "III/IV")) %>%
    filter(!is.na(FIGO_stage)) %>%
    ggplot(aes(x = FIGO_stage, y = mean_proba)) +
    geom_violin(aes(fill = FIGO_stage, color = FIGO_stage ), alpha = 0.7) +
    geom_point(aes(color = FIGO_stage )) +
    theme_classic() +
    ggpubr::stat_compare_means(comparisons = list(c("I/II", "III/IV")), label = "p.format") +
    labs(x = "", y = "Mean Cancer Scores") +
    scale_fill_manual(values = c("#4682B4", "#DC143C")) +
    scale_color_manual(values = c("#4682B4", "#DC143C")) +
    theme(legend.position = "none") +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


p_age <- df %>%
    filter(i.group == "OV") %>%
    mutate(age = if_else(Age >= 50, ">=50 yrs", "<50 yrs")) %>%
    group_by(sampleID, age) %>%
    summarise(mean_proba = mean(proba_1)) %>%
    ggplot(aes(x = age, y = mean_proba)) +
    geom_violin(aes(fill = age, color = age), alpha = 0.7) +
    geom_point(aes(color = age)) +
    theme_classic() +
    ggpubr::stat_compare_means(comparisons = list(c(">=50 yrs", "<50 yrs")),
                               label = "p.format") +
    labs(x = "", y = "Mean Cancer Scores") +
    scale_fill_manual(values = c("#4682B4", "#DC143C")) +
    scale_color_manual(values = c("#4682B4", "#DC143C")) +
    theme(legend.position = "none") 


df %>%
    filter(i.group == "OV") %>%
    mutate(age = if_else(Age >= 50, ">=50 yrs", "<50 yrs")) %>%
    group_by(sampleID, age) %>%
    summarise(mean_proba = mean(proba_1)) %>%
    ungroup() %>%
    group_by(age) %>%
    summarise(median = mean(mean_proba),
              sd = sd(mean_proba))


library(patchwork)
p <- p_age + p_figo + p_tfx

ggsave(filename = file.path(wrkdir,"./Fig3E_all_cancer_scores.png"), 
       p, width = 12, height = 5)

