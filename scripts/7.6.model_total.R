rm(list = ls())
library(metafor)
library(dplyr)
library(ggplot2)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("output/dataframes/effect_size.csv")
df <- df %>% filter(vi > 0)
df$replicate_ID[is.na(df$replicate_ID)] <- "0"
df$taxon_scientific[is.na(df$taxon_scientific)] <- "NA"
df$land_use_context_ur <- relevel(factor(df$land_use_context_ur), ref = "rural")
df$bd_metric_type <- as.factor(df$bd_metric_type)
df$land_use_context_details <- df$land_use_context_details/100
df$snh500 <- df$snh500/100
df$snh250 <- df$snh250/100
df$snh_CLC <- df$snh_CLC/100

df$cluster_ID <- interaction(
  df$study_ID,
  df$intervention_treatment_ID,
  df$control_treatment_ID,
  df$replicate_ID,
  df$taxon_other,                                                                                                                                 
  df$taxon_scientific,
  df$bd_metric_type,
  drop = TRUE
)  

df_abund <- df %>% filter(bd_metric_type == "abundance")
df_rich <- df %>% filter(bd_metric_type == "richness")
df_context <- df[!is.na(df$land_use_context_details), ]
df_context_CLC <- df_context[!is.na(df_context$snh_CLC), ]

df_size <- df[!is.na(df$intervention_size_sqm), ]
df_context_size <- df_context[!is.na(df_context$intervention_size_sqm), ]


##overall
overall <- rma.mv(yi, vi, random = ~1 | study_ID/cluster_ID, data = df, test = "t", method = "REML")
summary(overall)
overall_abund <- rma.mv(yi, vi, random = ~1 | study_ID/cluster_ID, data = df_abund, test = "t",method = "REML")
overall_rich <- rma.mv(yi, vi, random = ~1 | study_ID/cluster_ID, data = df_rich, test = "t",method = "REML")
summary(overall_abund)
summary(overall_rich)

#one-sided log-likelihood-ratio tests
overall_ml <- rma.mv(yi, vi, random = ~1 | study_ID/cluster_ID, data = df, test = "t", method = "ML")
mod_no_cluster <- rma.mv(yi, vi, random = ~1 | study_ID, data = df, test = "t", method = "ML")
anova(mod_no_cluster, overall_ml)
mod_no_study <- rma.mv(yi, vi, random = ~1 | cluster_ID, data = df, method = "ML") 
anova(mod_no_study, overall_ml)


##moderator_urban_rural
urban_rural <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_ur,
  random = ~1 | study_ID/cluster_ID,
  data = df,
  test = "t",
  method = "REML"
)
summary(urban_rural)

##moderator_urban_rural*metric_type
urban_rural_1 <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_ur*bd_metric_type,
  random = ~1 | study_ID/cluster_ID,
  data = df,
  test = "t",
  method = "REML"
)
summary(urban_rural_1)

##moderator_semi1000
total_semi1000 <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi1000)
AIC.rma(total_semi1000)

#moderator_semi1000_log
total_semi1000_log <- rma.mv(
  yi, vi,
  mods = ~ log(land_use_context_details),
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi1000_log)
AIC.rma(total_semi1000_log)

##moderator_semi1000_sq
df_context$snh1000_c <- scale(df_context$land_use_context_details, center = TRUE, scale = FALSE)
df_context$snh1000_c2 <- df_context$snh1000_c ^2
total_semi1000_sq <- rma.mv(
  yi, vi,
  mods = ~ snh1000_c + snh1000_c2,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
mean(df_context$land_use_context_details)
summary(total_semi1000_sq)
AIC.rma(total_semi1000_sq)

b1 <- coef(total_semi1000_sq)["snh1000_c"]
b2 <- coef(total_semi1000_sq)["snh1000_c2"]
vertex_c <- -b1 / (2 * b2)

snh_mean <- attr(df_context$snh1000_c, "scaled:center")  
vertex_original <- vertex_c + snh_mean

##moderator_snh1000*metric_type
total_semi1000_1 <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details*bd_metric_type,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi1000_1)

##moderator_size
total_size <- rma.mv(
  yi, vi,
  mods = ~ intervention_size_sqm,
  random = ~1 | study_ID/cluster_ID,
  data = df_size,
  test = "t",
  method = "REML"
)
summary(total_size)
AIC.rma(total_size)

#moderator_size_log
total_size_log <- rma.mv(
  yi, vi,
  mods = log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = df_size,
  test = "t",
  method = "REML"
)
AIC.rma(total_size_log)
summary(total_size_log)

#moderator_size_segmented500
log_threshold <- log(500)
df_size$log_size <- log(df_size$intervention_size_sqm)
df_size$log_below <- pmin(df_size$log_size, log_threshold)
df_size$log_above <- pmax(df_size$log_size - log_threshold, 0)

total_size_segmented_log <- rma.mv(
  yi, vi,
  mods = ~ log_below + log_above,
  random = ~1 | study_ID/cluster_ID,
  data = df_size,
  test = "t",
  method = "REML"
)
AIC.rma(total_size_segmented_log)
summary(total_size_segmented_log)

##moderator_snh+size
total_semi_size <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details + log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = df_context_size,
  test = "t",
  method = "REML"
)
summary(total_semi_size)

##moderator_semi1000*size
total_snh_size <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details*log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = df_context_size,
  test = "t",
  method = "REML"
)
summary(total_snh_size)

##moderator_ur+size
total_ur_size <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_ur+log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = df_size,
  test = "t",
  method = "REML"
)
summary(total_ur_size)


##sensitivity
##sensitivity_semi500
total_semi500 <- rma.mv(
  yi, vi,
  mods = ~ snh500,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi500)

##sensitivity_semi250
total_semi250 <- rma.mv(
  yi, vi,
  mods = ~ snh250,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi250)

##sensitivity_semi_euCLC
total_semi_CLC <- rma.mv(
  yi, vi,
  mods = ~ snh_CLC,
  random = ~1 | study_ID/cluster_ID,
  data = df_context_CLC,
  test = "t",
  method = "REML"
)
summary(total_semi_CLC)

##sensitivity_semi1000_euESA
total_semi1000_eu <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details,
  random = ~1 | study_ID/cluster_ID,
  data = df_context_CLC,
  test = "t",
  method = "REML"
)
summary(total_semi1000_eu)

##sensitivity_semi1000CV
total_semi1000CV <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details,
  random = ~1 | study_ID/cluster_ID,
  data = subset(df_context, context_cv <= 30),
  test = "t",
  method = "REML"
)
summary(total_semi1000CV)

#sensitivity_size_CV
total_size_logCV <- rma.mv(
  yi, vi,
  mods = log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = subset(df_size,size_cv <= 30),
  test = "t",
  method = "REML"
)
summary(total_size_logCV)

##publication bias
##multi-moderator multilevel meta-regression model
df$n_eff <- (df$i_sample_size * df$c_sample_size) / (df$i_sample_size + df$c_sample_size)
mod_sample_size <- rma.mv(yi, V = vi,
                          mods = ~ n_eff,
                          random = ~ 1 | study_ID/cluster_ID,
                          data = df, 
                          test = "t", 
                          method = "REML")

summary(mod_sample_size)

##Funnel plot
library(viridisLite)
blups <- ranef(overall)  
pos <- match(df$study_ID, unique(df$study_ID)) 
vi_study <- (blups$study_ID[[2]][pos])^2  
vi_residm <- df$vi  
vi_residc1 <- vi_residm + vi_study 
residc1 <- residuals(overall, type = "rstandard", level = 1)
par(
  cex.lab = 1.3,          
  cex.axis = 1,         
  cex.main = 1.3,        
  col.lab = "grey30",     
  col.axis = "grey30",    
  col.main = "grey30"    
)

funnel(
  residc1,                             
  vi = vi_residc1,                     
  yaxis = "sei",                      
  col = "grey30",                    
  xlim = c(-5, 5),                    
  ylab = "Standard error (SE)",       
  xlab = "Residuals conditional 1 (LRR)"  
)


## plots
##fig_2
library(ggtext)
res <- data.frame(
  label = c("Abundance (168/16)", "Richness (86/11)", "Overall Effect Size (254/18)"),
  est = c(coef(overall_abund), coef(overall_rich), coef(overall)),
  se = c(overall_abund$se, overall_rich$se, overall$se)
)
res$ci.lb <- res$est - 1.96 * res$se
res$ci.ub <- res$est + 1.96 * res$se
res$label_fmt <- sub(
  "^(.*?)\\s*(\\(.*\\))$",
  "<span style='font-weight:bold; color:black;'>\\1</span> \\2",
  res$label
)
res$label_fmt <- factor(label_split, levels = rev(label_split))  

overall_fig <- ggplot(res, aes(x = est, y = label_fmt)) +
  geom_point(size = 2, shape = 15,color = "grey30") +  
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = 0.15,,color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  xlab("LRR (± 95% CI)") +
  ylab(NULL) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),              
    axis.text.y = ggtext::element_markdown(size = 10), 
    axis.text.x = element_text(size = 10),     
    axis.title.x = element_text(size = 10, margin = margin(t = 15),,color = "grey30"),     
    plot.margin = margin(10, 30, 10, 10),       
    axis.ticks.y = element_blank()
  )

overall_fig
ggsave(filename = "output/figures/overall effect.png",plot = overall_fig,width = 5,height = 1.5,units = "in",dpi = 300,bg = "white")

##fig_3
new_data <- data.frame(land_use_context_details = seq(min(df_context$land_use_context_details),
                                                      max(df_context$land_use_context_details),
                                                      length.out = 100))
preds <- predict(total_semi1000, newmods = new_data$land_use_context_details, 
                 addx = TRUE, digits = 3)
plot_data <- data.frame(
  land_use_context_details = new_data$land_use_context_details,
  pred = preds$pred,
  ci.lb = preds$ci.lb,
  ci.ub = preds$ci.ub
)
snh_fig <- ggplot() +
  geom_point(data = df_context,
             aes(x = land_use_context_details, y = yi, size = 1/vi),
             color = "grey30", alpha = 0.5) + 
  geom_line(data = plot_data,
            aes(x = land_use_context_details, y = pred),
            color = "cadetblue", size = 1) +
  geom_ribbon(data = plot_data,
              aes(x = land_use_context_details, ymin = ci.lb, ymax = ci.ub),
              fill = "cadetblue", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Proportion of Semi-Natural Habitat",
       y = "Predicted Effect Size (LRR)",
       size = "Weight (1/vi)") +
  scale_size_continuous(range = c(1, 4)) +
  theme(
    panel.background = element_blank(),                                 
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3), 
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5), 
    legend.position = "bottom",
    legend.title = element_text(size = 10, color = "grey30"),
    legend.text  = element_text(size = 10, color = "grey30"),
    axis.title = element_text(size = 12, color = "grey30"),
    axis.text = element_text(size = 10, color = "grey30")
  )
snh_fig
ggsave(filename = "output/figures/snh.png",plot = snh_fig,width = 5,height = 4,units = "in",dpi = 300,bg = "white")

#fig_A1
group_n <- df %>%
  group_by(land_use_context_ur, bd_metric_type) %>%
  summarise(n = n(), .groups = "drop")

new_data <- data.frame(
  land_use_context_ur = rep(c("rural", "urban"), each = 2),
  bd_metric_type = rep(c("abundance", "richness"), 2)
)
mod_matrix <- model.matrix(~ land_use_context_ur * bd_metric_type, data = new_data)[, -1]
pred <- predict(urban_rural_1, newmods = mod_matrix)
plot_data <- new_data %>%
  left_join(group_n, by = c("land_use_context_ur", "bd_metric_type")) %>%
  mutate(
    Estimate = pred$pred,
    CI_lower = pred$ci.lb,
    CI_upper = pred$ci.ub,
    Position =  c(1.5, 2, 3, 3.5),              
    Shape = rep(c("abundance", "richness"), 2)   
  )

shape_map <- c("abundance" = 16, "richness" = 21)
ur_fig <- ggplot(plot_data, aes(x = Position, y = Estimate)) +
  geom_segment(aes(xend = Position, y = CI_lower, yend = CI_upper), color = "grey30", linewidth = 0.3) +
  
  geom_point(aes(shape = Shape), size = 2, fill = "white", stroke = 0.8, color = "grey30") +
  scale_shape_manual(name = "",values = shape_map) +
  geom_text(aes(label = n, y = CI_upper + 0.3), size = 3.2, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  scale_x_continuous(
    breaks =  c(1.75, 3.25),
    labels = c("Rural", "Urban")
  ) +
  coord_cartesian(xlim = c(1, 3.8))+
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 11, vjust = 1,color = "grey30"),
    axis.line.y = element_line(color = "grey30", linewidth = 0.3),
    axis.title.y = element_text(size = 11, color = "grey30", margin = margin(r = 10)),
    legend.position = "bottom", 
    legend.box.margin = margin(t = -25),
    legend.title = element_text(size = 11,color = "grey30"),
    legend.text = element_text(size = 11,color = "grey30")
  ) +
  labs(
    x = "",
    y = "LRR (± 95% CI)"
  )
ur_fig
ggsave(filename = "output/figures/ur.png",plot = ur_fig,width = 5,height = 3,units = "in",dpi = 300,bg = "white")

#fig_A2
new_data_1 <- data.frame(
  land_use_context_details = seq(
    min(df_context$land_use_context_details, na.rm = TRUE),
    max(df_context$land_use_context_details, na.rm = TRUE),
    length.out = 100
  )
)
new_data_1$snh1000_c  <- new_data_1$land_use_context_details - mean(df_context$land_use_context_details, na.rm = TRUE)
new_data_1$snh1000_c2 <- new_data_1$snh1000_c^2
pred <- predict(total_semi1000_sq, newmods = as.matrix(new_data_1[, c("snh1000_c", "snh1000_c2")]))
plot_data <- cbind(new_data_1, pred)

snh_sq_fig <- ggplot(plot_data, aes(x = land_use_context_details, y = pred)) +
  geom_point(data = df_context,
             aes(x = land_use_context_details, y = yi, size = 1/vi),
             color = "grey30", alpha = 0.5) +
  scale_size_continuous(range = c(1, 4), name = "Weight (1/vi)") +
  geom_line(color = "cadetblue", size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "cadetblue", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x = "Proportion of Semi-Natural Habitat",
    y = "Predicted Effect Size (LRR)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    axis.title = element_text(size = 12, color = "grey30"),
    axis.text = element_text(size = 10, color = "grey30"),
    legend.text = element_text(size = 10, color = "grey30"),
    legend.title = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )
snh_sq_fig
ggsave(filename = "output/figures/snh_sq.png",plot = snh_sq_fig,width = 5,height = 4,units = "in",dpi = 300,bg = "white")


##residual plots
#overall
res_overall <- residuals(overall)
qqnorm(res_overall)
qqline(res_overall, col = "red")
fitted_vals_overall <- fitted(overall)
plot(fitted_vals_overall, res_overall,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
#ur
res_ur <- residuals(urban_rural)
qqnorm(res_ur)
qqline(res_ur , col = "red")
fitted_vals_ur <- fitted(urban_rural)
plot(fitted_vals_ur, res_ur,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

#snh
res_SNH <- residuals(total_semi1000)
qqnorm(res_SNH )
qqline(res_SNH , col = "red")
fitted_vals <- fitted(total_semi1000)
plot(fitted_vals, res_SNH,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

#snh_sq 
res_SNH_sq <- residuals(total_semi1000_sq)
qqnorm(res_SNH_sq)
qqline(res_SNH_sq, col = "red")
fitted_vals_sq <- fitted(total_semi1000_sq)
plot(fitted_vals_sq,res_SNH_sq,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

#size
res_size <- residuals(total_size )
qqnorm(res_size)
qqline(res_size, col = "red")
fitted_vals_size <- fitted(total_size )
plot(fitted_vals_size, res_size,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

#size_log
res_size_log <- residuals(total_size_log)
hist(res_size_log, main = "Histogram of residuals", xlab = "Residuals")
hist(df_size$yi)
qqnorm(res_size_log)
qqline(res_size_log, col = "red")
fitted_vals_size_log <- fitted(total_size_log)
plot(fitted_vals_size_log,res_size_log,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")


