## librries
library(dplyr)
library(ggplot2)

## load data
dir <- "C:/PredictioR/result/FL/integrateAI_04122024"

res.pan <- read.csv(file.path(dir, "pan_cancer.csv"))
res.percan <- read.csv(file.path(dir, "per_cancer.csv"))
res.pertreatment <- read.csv(file.path(dir, "per_treatment.csv"))

##################################################################################
################################### per cancer ###################################
##################################################################################

## Per cancer: visualization plot
group <- unique(res.percan$Cancer_type)

res <- lapply(1:length(group), function(i){

  df <- res.percan[res.percan$Cancer_type == group[i], ]
  data.frame( cancer_type = df$Cancer_type[1],
              gene = df$Gene,
              method = c(rep("FL", nrow(df)),
                         rep("FE", nrow(df)),
                         rep("RE", nrow(df))),
              coef = c( df$coef.FL,
                        df$Coef.fixed,
                        df$Coef),
              se =  c( df$se.FL,
                       df$SE.fixed,
                       df$SE),
              pval =  c( df$Pval.FL,
                       df$Pval.fixed,
                       df$Pval),
              I2 =  c( rep(NA, nrow(df)),
                       df$I2,
                       df$I2),
              Q_pval =  c( rep(NA, nrow(df)),
                       df$Q_Pval,
                       df$Q_Pval))


})

res <- do.call(rbind, res)

## compare distributions using estimated integrated coef (FL vs FE)
res.percan$I2_group <-  ifelse(res.percan$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "scatterplot_percan_coef_FL_RE.jpeg"),
     width = 800, height = 600, res = 150)

ggplot(data = res.percan) +
  geom_point(mapping = aes(x = coef.FL, y = Coef, colour = I2_group)) +
  facet_grid(~ Cancer_type) +
  labs(x="integrated coef (FL)", y="integrated coef (RE)" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## scatterplot using estimated integrated coef
jpeg(file=file.path(dir, "scatterplot_percan_coef.jpeg"),
     width = 800, height = 600, res = 150)

res %>%
  ggplot(aes(method, coef, color = method, fill = method)) +
  geom_violin() +
  facet_grid(~ cancer_type) +
  labs(x=NULL, y="integrated coef" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## compare distributions using estimated integrated pvalue
jpeg(file=file.path(dir, "violinplot_percan_pval.jpeg"),
     width = 800, height = 600, res = 150)

res %>%
  ggplot(aes(method, pval, color = method, fill = method)) +
  geom_violin() +
  facet_grid(~ cancer_type) +
  labs(x=NULL, y="integrated pval" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))

dev.off()

## compare heterogeneity

rem.FL <- res[!is.na(res$I2),]
rem.FL$I2_group <- ifelse(rem.FL$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "violinplot_percan_I2.jpeg"),
     width = 800, height = 600, res = 150)

rem.FL %>%
  ggplot(aes(I2_group, I2, color = I2_group, fill = I2_group)) +
  geom_violin() +
  facet_grid(~ cancer_type) +
  labs(x=NULL, y="heterogeneity I2" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## correlation test
res <- lapply(1:length(group), function(i){

  df <- res.percan[res.percan$Cancer_type == group[i], ]
  fit.fl.fe <- cor.test(df$coef.FL, df$Coef.fixed)
  fit.fl.re <- cor.test(df$coef.FL, df$Coef)
  res <- data.frame(cancer_type = group[i],
                    r_fl_fe = fit.fl.fe$estimate,
                    pval_fl_fe = fit.fl.fe$p.value,
                    r_fl_re = fit.fl.re$estimate,
                    pval_fl_re = fit.fl.re$p.value,
                    I2_lower = nrow(df[df$I2 < 0.50, ]),
                    I2_upper = nrow(df[df$I2 >= 0.50, ]))


})

res <- do.call(rbind, res)
write.csv(res, file.path(dir, "cor_percan.csv" ), row.names = FALSE)


##################################################################################
################################### per treatment ################################
##################################################################################

## Per treatment: visualization plot
group <- unique(res.pertreatment$Treatment)

res <- lapply(1:length(group), function(i){

  df <- res.pertreatment[res.pertreatment$Treatment == group[i], ]
  data.frame( treatment_type = df$Treatment[1],
              method = c(rep("FL", nrow(df)),
                         rep("FE", nrow(df)),
                         rep("RE", nrow(df))),
              coef = c( df$coef.FL,
                        df$Coef.fixed,
                        df$Coef),
              se =  c( df$se.FL,
                       df$SE.fixed,
                       df$SE),
              pval =  c( df$Pval.FL,
                         df$Pval.fixed,
                         df$Pval),
              I2 =  c( rep(NA, nrow(df)),
                       df$I2,
                       df$I2))


})

res <- do.call(rbind, res)

## compare distributions using estimated integrated coef (FL vs FE)
res.pertreatment$I2_group <-  ifelse(res.pertreatment$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "scatterplot_pertreatment_coef_FL_FE.jpeg"),
     width = 800, height = 600, res = 150)

ggplot(data = res.pertreatment) +
  geom_point(mapping = aes(x = coef.FL, y = Coef.fixed, colour = I2_group)) +
  facet_grid(~ Treatment) +
  labs(x="integrated coef (FL)", y="integrated coef (FE)" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8, face="bold"))


dev.off()
## compare distributions
jpeg(file=file.path(dir, "violinplot_pertreatment_coef.jpeg"),
     width = 600, height = 600, res = 150)

res %>%
  ggplot(aes(method, coef, color = method, fill = method)) +
  geom_violin() +
  facet_grid(~ treatment_type) +
  labs(x=NULL, y="integrated coef" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## compare heterogeneity

rem.FL <- res[!is.na(res$I2),]
rem.FL$I2_group <- ifelse(rem.FL$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "violinplot_pertreatment_I2.jpeg"),
     width = 800, height = 600, res = 150)

rem.FL %>%
  ggplot(aes(I2_group, I2, color = I2_group, fill = I2_group)) +
  geom_violin() +
  facet_grid(~ treatment_type) +
  labs(x=NULL, y="heterogeneity I2" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## correlation test
res <- lapply(1:length(group), function(i){

  df <- res.pertreatment[res.pertreatment$Treatment == group[i], ]
  fit.fl.fe <- cor.test(df$coef.FL, df$Coef.fixed)
  fit.fl.re <- cor.test(df$coef.FL, df$Coef)
  res <- data.frame(treatment_type = group[i],
                    r_fl_fe = fit.fl.fe$estimate,
                    pval_fl_fe = fit.fl.fe$p.value,
                    r_fl_re = fit.fl.re$estimate,
                    pval_fl_re = fit.fl.re$p.value,
                    I2_lower = nrow(df[df$I2 < 0.50, ]),
                    I2_upper = nrow(df[df$I2 >= 0.50, ]))

})

res <- do.call(rbind, res)
write.csv(res, file.path(dir, "cor_pertreatment.csv" ), row.names = FALSE)


##################################################################################
################################### pan cancer ################################
##################################################################################

## Pan cancer: visualization plot
group <- c("FL", "FE", "RE")
df <- res.pan
res <- lapply(1:length(group), function(i){

  data.frame( method = c(rep("FL", nrow(df)),
                         rep("FE", nrow(df)),
                         rep("RE", nrow(df))),
              coef = c( df$coef.FL,
                        df$Coef.fixed,
                        df$Coef),
              se =  c( df$se.FL,
                       df$SE.fixed,
                       df$SE),
              pval =  c( df$Pval.FL,
                         df$Pval.fixed,
                         df$Pval),
              I2 =  c( rep(NA, nrow(df)),
                       df$I2,
                       df$I2))


})

res <- do.call(rbind, res)

## compare distributions
jpeg(file=file.path(dir, "violinplot_pan_coef.jpeg"),
     width = 600, height = 600, res = 150)

res %>%
  ggplot(aes(method, coef, color = method, fill = method)) +
  geom_violin() +
  #facet_grid(~ treatment_type) +
  labs(x=NULL, y="integrated coef" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## compare heterogeneity

rem.FL <- res[!is.na(res$I2),]
rem.FL$I2_group <- ifelse(rem.FL$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "violinplot_pan_I2.jpeg"),
     width = 800, height = 600, res = 150)

rem.FL %>%
  ggplot(aes(I2_group, I2, color = I2_group, fill = I2_group)) +
  geom_violin() +
  #facet_grid(~ treatment_type) +
  labs(x=NULL, y="heterogeneity I2" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="none",
    legend.text = element_text(size = 8, face="bold"))


dev.off()

## correlation test

df <- res.pan[res.pan$Pval < 0.05, ]
fit.fl.fe <- cor.test(df$coef.FL, df$Coef.fixed)
fit.fl.re <- cor.test(df$coef.FL, df$Coef)
res <- data.frame(r_fl_fe = fit.fl.fe$estimate,
           pval_fl_fe = fit.fl.fe$p.value,
           r_fl_re = fit.fl.re$estimate,
           pval_fl_re = fit.fl.re$p.value,
           I2_lower = nrow(df[df$I2 < 0.50, ]),
           I2_upper = nrow(df[df$I2 >= 0.50, ]))

write.csv(res, file.path(dir, "cor_pan.csv" ), row.names = FALSE)


## compare distributions using estimated integrated coef (FL vs FE)
res.pan$I2_group <-  ifelse(res.pan$I2 > 0.50, "I2 > 50%", "I2 <= 50%")

jpeg(file=file.path(dir, "scatterplot_pan_coef_FL_FE.jpeg"),
     width = 800, height = 600, res = 150)

ggplot(data = res.pan) +
  geom_point(mapping = aes(x = coef.FL, y = Coef.fixed, colour = I2_group)) +
  #facet_grid(~ Treatment) +
  labs(x="integrated coef (FL)", y="integrated coef (FE)" ) +
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8, face="bold"))


dev.off()
