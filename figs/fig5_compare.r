library(tidyverse)
library(readxl)
library(data.table)
library(pROC)


df <- read_excel("./ROMA_vs_frag.csv")

model_combined <- roc(df$group, df$mean_comb_proba)
model_frag <- roc(df$group, df$mean_proba)
model_ca125 <- roc(df$group, df$CA125)
model_ca125_pred <- roc(df$group, df$CA125_pred)
model_he4 <- roc(df$group, df$HE4)
model_he4_pred <- roc(df$group, df$HE4_pred)
model_roma <- roc(df$group, df$ROMA)
model_roma_pred <- roc(df$group, df$ROMA_pred)


# auc_combined <- round(auc(df$group, df$mean_comb_proba)[1], 3)
auc_frag <- round(auc(df$group, df$mean_proba)[1], 3)
auc_ca125 <- round(auc(df$group, df$CA125)[1], 3)
auc_ca125_pred <- round(auc(df$group, df$CA125_pred)[1], 3)
auc_he4 <- round(auc(df$group, df$HE4)[1], 3)
auc_he4_pred <- round(auc(df$group, df$HE4_pred)[1], 3)
auc_roma <- round(auc(df$group, df$ROMA)[1], 3)
auc_roma_pred <- round(auc(df$group, df$ROMA_pred)[1], 3)

# ci_combined <- round(ci.auc(df$group, df$mean_comb_proba), 3)
ci_frag <- round(ci.auc(df$group, df$mean_proba), 3)
ci_ca125 <- round(ci.auc(df$group, df$CA125), 3)
ci_ca125_pred <- round(ci.auc(df$group, df$CA125_pred), 3)
ci_he4 <- round(ci.auc(df$group, df$HE4), 3)
ci_he4_pred <- round(ci.auc(df$group, df$HE4_pred), 3)
ci_roma <- round(ci.auc(df$group, df$ROMA), 3)
ci_roma_pred <- round(ci.auc(df$group, df$ROMA_pred), 3)


pdf("./Fig5_compare.pdf", 
    height = 7, width = 7)

plot(model_frag, col = "black", 
     xlab = "False Positive Rate", 
     ylab = "True Positive Rate",
     asp = NA,
     grid = TRUE,
     legacy.axes = TRUE)

plot(model_ca125, col = "#FF7F50", lty = "dotted", add = T)
plot(model_ca125_pred, col = "#FF7F50", add = T)
plot(model_he4, col = "#20B2AA", lty = "dotted", add = T)
plot(model_he4_pred, col = "#20B2AA", add = T)
plot(model_roma, col = "#9ACD32", lty = "dotted", add = T)
plot(model_roma_pred, col = "#9ACD32", add = T)


legend("bottomright",
       legend = c(paste0("Fragmentomics: ", auc_frag, ", 95%CI:", ci_frag[1], "-", ci_frag[3]), 
                  paste0("CA125(value): ", auc_ca125, ", 95%CI:", ci_ca125[1], "-", ci_ca125[3]),
                  paste0("CA125: ", auc_ca125_pred, ", 95%CI:", ci_ca125_pred[1], "-", ci_ca125_pred[3]),
                  paste0("HE4(value): ", auc_he4, ", 95%CI:", ci_he4[1], "-", ci_he4[3]),
                  paste0("HE4: ", auc_he4_pred, ", 95%CI:", ci_he4[1], "-", ci_he4_pred[3]),
                  paste0("ROMA index(value): ", auc_roma, ", 95%CI:", ci_roma[1], "-", ci_roma[3]),
                  paste0("ROMA index: ", auc_roma_pred, ", 95%CI:", ci_roma_pred[1], "-", ci_roma_pred[3])),
       col = c("black", "#FF7F50", "#FF7F50", "#20B2AA", "#20B2AA", "#9ACD32", "#9ACD32"),
       lty = c("solid", "dotted", "solid", "dotted", "solid", "dotted", "solid"),
       lwd=2, cex = 0.8, xpd = TRUE)

dev.off()