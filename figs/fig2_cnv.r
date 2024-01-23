library(data.table)
library(tidyverse)
library(biomaRt)
library(pacman)



# fig 2a.
# ---------------------------------------------------------
# get ensembl
ensembl <- useMart("ensembl", 
                   dataset = "hsapiens_gene_ensembl")

ensembl@host # "https://www.ensembl.org:443/biomart/martservice"


# annotate genes
gene_cnv <- fread("./TCGA.OV.sampleMap_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz")
all_genes <- unique(gene_cnv$`Gene Symbol`)

gene_info <- getBM(attributes = c("chromosome_name", 
                                  "start_position", 
                                  "end_position", 
                                  "strand",
                                  "hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = all_genes,
                   mart = ensembl)

merged <- merge.data.table(gene_cnv, gene_info,
                           by.x = "Gene Symbol", 
                           by.y = "hgnc_symbol")


med_cohort <- merged %>%
    filter(chromosome_name %in% seq(1, 22)) %>%
    dplyr::select(chromosome_name, start_position, starts_with("TCGA")) %>%
    rowwise() %>%
    mutate(med = median(c_across(starts_with("TCGA"))),
           std = sd(c_across(starts_with("TCGA")))) %>% 
    ungroup() %>%
    dplyr::select(chromosome_name, start_position, med, std)


med_cohort %>%
    mutate(chromosome_name = factor(chromosome_name, level = seq(1,22)),
           cnv_type = if_else(med <= 0 , "del", "amp")) %>%
    ggplot(aes(x = start_position, y = med)) +
    geom_line(aes(color = cnv_type)) +
    geom_area(aes(fill = cnv_type), 
              alpha = 0.6) +
    facet_grid(~chromosome_name,
               space = "free_x",
               scales = "free_x") +
    labs(x = "", y = "Median LogR", 
         title = "TCGA OV cohort (n = 579)") +
    geom_hline(yintercept = 0) + 
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0, "pt"),
          legend.position="none",
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#E3170D", "#191970")) +
    scale_color_manual(values = c("#E3170D", "#191970"))

ggsave("./TCGA_OV_CNV.pdf", width = 9, height = 4)



# fig 2b.
# ---------------------------------------------------------
library(ComplexHeatmap)
library(readxl)

logr <- read.csv("../../results/features/CNV.csv", 
                 check.names = F)
samples <- read_excel("./Supplementary_Tables.xlsx")
logr <- subset(logr, select = c("chr", "start", "end", samples$sampleID))

ov <- samples %>% 
    filter(group == "OV") %>%
    left_join(tfx, by = "sampleID") %>% 
    arrange(-tfx) %>%
    pull(sampleID)

health <- samples %>% 
    filter(group == "Healths") %>% 
    left_join(tfx, by = "sampleID") %>%
    arrange(-tfx) %>%
    pull(sampleID)

benign <- samples %>% 
    filter(group == "Benigns") %>% 
    left_join(tfx, by = "sampleID") %>%
    arrange(-tfx) %>%
    pull(sampleID)


logr_mat <- logr %>%
    mutate(label = paste0(chr, "-", start)) %>%
    select(-c(1:3)) %>%
    select(label, everything()) %>%
    select(label, ov, benign, health) %>%
    rowwise() %>%
    mutate(mean_logr = mean(c_across(starts_with("MRD")))) %>%
    filter(mean_logr > -1) %>%
    column_to_rownames("label") %>%
    select(-mean_logr,) %>%
    as.matrix()  %>%
    t()

groups <- c(rep("OV", length(ov)), 
            rep("Benign", length(benign)),
            rep("Healthy", length(health)))

chroms <- logr %>%
    rowwise() %>%
    mutate(mean_logr = mean(c_across(starts_with("MRD")))) %>%
    filter(mean_logr > -1) %>%
    pull(chr)


chrom_anno <- HeatmapAnnotation(chrom = anno_block(gp = gpar(fill = c("#D3D3D3", "#F5F5F5"),
                                                             lineheight = 3),
                                                   labels = seq(1,22),
                                                   labels_gp = gpar(fontsize = 8)))


# svg(filename = "CNV_heatmap.svg", width = 18, height = 8)
pdf(width = 9.6, height = 4.5)
ComplexHeatmap::Heatmap(logr_mat, 
                        cluster_columns = F,
                        cluster_rows = F,
                        show_row_names = F,
                        show_column_names = F,
                        row_split = groups,
                        column_split = chroms,
                        top_annotation = chrom_anno,
                        cluster_row_slices = F,
                        cluster_column_slices = F,
                        column_title = NULL,
                        heatmap_legend_param = list(title = "logR ratio"))

dev.off()


# fig 2c.
# -------------------------------------------
# correlation between tcga logR and pumc logR
tfx <- fread("../../results/features/ichorCNA_tfx.tsv")
colnames(tfx) <- c("sampleID", "tfx", "ploidy")
tfx[, sampleID:=gsub("_tumor", "", sampleID)]


df <- tfx %>%
    right_join(samples, by = "sampleID") 

high_tfx <- df %>%
    filter(tfx >= 0.03) %>%
    pull(sampleID)

pumc_logR <- logr %>%
    select(chr, start, end, high_tfx) %>%
    rowwise() %>%
    mutate(med_logR = mean(c_across(high_tfx))) %>%
    select(chr, start, med_logR)


tcga_logR <- fread("../..//manuscript/TCGA_OV_CNV_plotting/median_logR.csv")

tcga_logR <- tcga_logR %>%
    mutate(start = round(start_position / 1000000) * 1000000 + 1) %>%
    select(chromosome_name, start, med) %>%
    rename("chr" = "chromosome_name")

join <- pumc_logR %>%
    filter(med_logR > -0.5) %>%
    inner_join(tcga_logR, by = c("chr", "start"))

cor.test(join$med, join$med_logR)


ggplot(join, aes(x = med, y = med_logR)) +
    geom_point() +
    geom_smooth(method = "lm") + 
    theme_bw() +
    labs(x = "median logR of TCGA OV", 
         y = "median logR of PUMC samples",
         title = "c") +
    annotate("text",
             label = "R = 0.65, p < 2.2e-16", 
             x = -0.4, y = 0.17)


ggsave("./Fig2C_cor.png", width = 5, height = 5)



# fig 2d.
# -------------------------------------------

df %>%
    dplyr::mutate(group = case_when(FIGO_stage == "I" ~ "I/II",
                             FIGO_stage == "II" ~ "I/II",
                             FIGO_stage == "III" ~ "III/IV",
                             FIGO_stage == "IV" ~ "III/IV",
                             Malignancy == "Benign" ~ "control",
                             Malignancy == "healthy" ~ "control")) %>%
    filter(!is.na(group)) %>%
    # group_by(group) %>%
    # summarise(median = mean(tfx))
    ggplot(aes(x = group, y = tfx)) +
    # geom_point() +
    geom_violin(aes(fill = group)) +
    ggpubr::stat_compare_means(comparisons = list(c("I/II", "control"),
                                                  c("III/IV", "control")),
                               method = "wilcox.test",
                               label = "p.signif") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", y = "Tumor fraction", title = "d") +
    scale_fill_manual(values = c("#6495ED", "#FA8072", "#DC143C")) +
    scale_y_continuous(labels = scales::percent)


ggsave("./Fig1D_tfx.png", width = 5, height = 5)

