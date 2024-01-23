library(tidyverse)
library(data.table)
library(readxl)

options(scipen=999)



# fig 3a.
# ---------------------------------------------------------

fsd <- fread("../../results/features/FSD.csv")
samples <- read_excel("./Supplementary_Tables.xlsx")


comb_fsd <- fsd %>%
    right_join(samples, by = "sampleID") %>%
    mutate(chrom = paste0(chrom, arm)) %>%
    # filter(chrom == "chr16p") %>%
    group_by(sampleID, start, group) %>%
    summarise(sum_count = sum(count)) %>%
    ungroup() %>%
    group_by(sampleID) %>%
    mutate(ratio = sum_count / sum(sum_count)) 


p_fsd <- comb_fsd %>%
    mutate(start = factor(start)) %>%
    mutate(groups = case_when(group == "OV" ~ "OV",
                              group == "Benigns" ~ "Benign",
                              group == "Healths" ~ "Healthy")) %>%
    mutate(groups = factor(groups, level = c("OV", "Benign", "Healthy"))) %>%
    ggplot() +
    geom_boxplot(aes(x = start, y = ratio, 
                     fill = groups, color = groups),
                 outlier.size = 1.2) +
    theme_bw() +
    scale_fill_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_color_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_x_discrete(breaks = factor(seq(50, 400, 50))) +
    labs(x = "", y = "Ratio") +
    theme(legend.position = c(0.87, 0.85),
          legend.title = element_blank())

ggsave("./fsd_dist.png", p_fsd, width = 10, height = 4)



# fig 3b.
# ---------------------------------------------------------

fsr <- fread("../../results/features//FSR.csv")


breaks <- fsr %>% 
    group_by(chr) %>% 
    filter(row_number() == ceiling(n()/2)) %>% 
    pull(bin)

labels <- fsr %>% 
    group_by(chr) %>% 
    filter(row_number() == ceiling(n()/2))  %>% 
    mutate(chr = gsub("chr", "", chr)) %>%
    pull(chr)


p_fsr <- fsr %>%
    right_join(samples, by = c("id" = "sampleID")) %>% 
    mutate(groups = case_when(group == "OV" ~ "OV",
                              group == "Benigns" ~ "Benign",
                              group == "Healths" ~ "Healthy")) %>%
    mutate(groups = factor(groups, level = c("OV", "Benign", "Healthy"))) %>%
    ggplot(aes(x = bin, y = ratio.centered, group = id)) +
    geom_line(aes(color = groups), size = .3, alpha = .3) +
    facet_wrap(~groups, dir = "v", 
               strip.position = "top") +
    theme_classic() +
    scale_fill_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_color_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    labs(x = "", y = "Centroid FSR") +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank()) 
    # scale_x_continuous(breaks = breaks, labels = labels) 

ggsave(filename = "./Fig2C_FSR.png",
       width = 10, height = 4)



# fig 3c.d.
# ---------------------------------------------------------

bpm <- fread("./BPM_unenrich.csv")
edm <- fread("./EDM_unenrich.csv")

bpm <- bpm %>% right_join(samples, by = "sampleID") 
edm <- edm %>% right_join(samples, by = "sampleID") 


# bpm
bpm_diff <- data.frame()
i <- 1
for (imotif in unique(bpm$motif)) {
    ov <- bpm %>% filter(motif == imotif, group == "OV") %>% pull(bpm_ratio)
    benign <- bpm %>% filter(motif == imotif, group == "Benigns") %>% pull(bpm_ratio)
    healths <- bpm %>% filter(motif == imotif, group == "Healths") %>% pull(bpm_ratio)
    ctrl <- bpm %>% filter(motif == imotif, group != "OV") %>% pull(bpm_ratio)
    
    p <- wilcox.test(ov, benign)$p.value
    bpm_diff[i, "motif"] <- imotif
    bpm_diff[i, "p_benign"] <- p
    p <- wilcox.test(ov, healths)$p.value
    bpm_diff[i, "p_healths"] <- p
    p <- wilcox.test(benign, healths)$p.value
    bpm_diff[i, "p_within_ctrl"] <- p
    
    bpm_diff[i, "med_ov"] <- median(ov)
    bpm_diff[i, "med_benign"] <- median(benign)
    bpm_diff[i, "med_health"] <- median(healths)
    i <- i + 1
}


bpm_diff <- bpm_diff %>% 
    as_tibble() %>%
    mutate_at(vars(starts_with("p_")), ~ p.adjust(.x, n = nrow(bpm_diff)))


p_bpm <- bpm_diff %>% 
    filter(p_healths < 0.05, p_benign < 0.05, p_within_ctrl > 0.05) %>%
    head(n = 6) %>%
    left_join(bpm, by = "motif") %>%
    mutate(group = factor(group, level = c("OV", "Benigns", "Healths"))) %>%
    ggplot(aes(x = group, y = bpm_ratio,
               fill = group, color = group)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~motif) +
    theme_bw() +
    scale_fill_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_color_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    ggpubr::stat_compare_means(comparisons = list(c("Benigns", "Healths"),
                                                  c("OV", "Benigns"),
                                                  c("OV", "Healths")),
                               label = "p.signif", vjust = 1.5) +
    labs(x = "", y = "BPM ratio") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("OV", "Benign", "Healthy"))

ggsave(filename = "../../manuscript/Fig2C_BPM.png", 
       p_bpm, width = 7, height = 7)


# edm
edm_diff <- data.frame()
i <- 1
for (imotif in unique(edm$motif)) {
    ov <- edm %>% filter(motif == imotif, group == "OV") %>% pull(edm_ratio)
    benign <- edm %>% filter(motif == imotif, group == "Benigns") %>% pull(edm_ratio)
    healths <- edm %>% filter(motif == imotif, group == "Healths") %>% pull(edm_ratio)
    ctrl <- edm %>% filter(motif == imotif, group != "OV") %>% pull(edm_ratio)
    
    p <- wilcox.test(ov, benign)$p.value
    edm_diff[i, "motif"] <- imotif
    edm_diff[i, "p_benign"] <- p
    p <- wilcox.test(ov, healths)$p.value
    edm_diff[i, "p_healths"] <- p
    p <- wilcox.test(benign, healths)$p.value
    edm_diff[i, "p_within_ctrl"] <- p
    
    edm_diff[i, "med_ov"] <- median(ov)
    edm_diff[i, "med_benign"] <- median(benign)
    edm_diff[i, "med_health"] <- median(healths)
    i <- i + 1
}


edm_diff <- edm_diff %>% 
    as_tibble() %>%
    mutate_at(vars(starts_with("p_")), ~ p.adjust(.x, n = nrow(edm_diff)))

p_edm <- edm_diff %>% 
    filter(p_healths < 0.05, p_benign < 0.05, p_within_ctrl > 0.05) %>% 
    head(n = 6) %>%
    left_join(edm, by = "motif") %>%
    mutate(group = factor(group, level = c("OV", "Benigns", "Healths"))) %>%
    ggplot(aes(x = group, y = edm_ratio,
               fill = group, color = group)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~motif) +
    theme_bw() +
    scale_fill_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    scale_color_manual(values = c("#DC143C", "#4682B4", "#87CEFA")) +
    ggpubr::stat_compare_means(comparisons = list(c("Benigns", "Healths"),
                                                  c("OV", "Benigns"),
                                                  c("OV", "Healths")),
                               label = "p.signif", vjust = 1.5) +
    labs(x = "", y = "EDM ratio") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("OV", "Benign", "Healthy"))


ggsave(filename = "../../manuscript/Fig2C_EDM.png", 
       p_edm,
       width = 7, height = 7)
