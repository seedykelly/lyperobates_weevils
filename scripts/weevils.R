# if git is ahead by X commits do this: git reset --soft HEAD~1 (8=# of commits)

## ---- analysis ----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(flextable)
library(officer)
library(kableExtra)
library(knitr)
library(officedown)
library(tidybayes)
library(gsg)
library(mgcv)
library(MASS)
library(tibble)
library(stringr)

# install.packages("~/Documents/Action/gsg_2.0.tar", repos = NULL, type = "source")

#### read in data file ####
morph_data <- read.csv(file="data/raw/weevils.csv", header=TRUE, sep=",", dec=".") %>%
  as_tibble() %>%
  mutate(fem=rowMeans(dplyr::select(., l_fem, r_fem), na.rm = TRUE)) %>%
  mutate(tib=rowMeans(dplyr::select(., l_tib, r_tib), na.rm = TRUE)) %>%
  mutate(total_body=tot_abdo+thorax) %>%
  mutate(total_leg=fem+tib) 

# =============================
# sexual size dimorphism
# =============================
# make sure sex is a factor
morph_data$sex <- factor(morph_data$sex)

morph_data %>%
  filter(sex=="m") %>%
  filter(total_leg<6)

summary_stats <- morph_data %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    
    mean_total_body = mean(total_body, na.rm = TRUE),
    se_total_body = sd(total_body, na.rm = TRUE) / sqrt(n),
    
    mean_total_leg = mean(total_leg, na.rm = TRUE),
    se_total_leg = sd(total_leg, na.rm = TRUE) / sqrt(n)
  )

summary_stats

# total body length
model_body <- lm(total_body ~ sex, data = morph_data)
summary(model_body)

# total leg length
model_leg <- lm(total_leg ~ sex, data = morph_data)
summary(model_leg)

manova_model <- manova(cbind(total_body, total_leg) ~ sex, data = morph_data)
summary(manova_model)

plot_data <- morph_data %>%
  dplyr::select(sex, total_body, total_leg) %>%
  pivot_longer(
    cols = c(total_body, total_leg),
    names_to = "trait",
    values_to = "value"
  )

figure_2 <- ggplot(plot_data, aes(x = trait, y = value, fill = sex)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.6,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.7
    ),
    size = 1.6,
    alpha = 0.7,
    color = "black"
  ) +
  scale_fill_grey(
    start = 0.8,
    end = 0.4,
    labels = c("f" = "Females", "m" = "Males")
  ) +
  scale_x_discrete(
    labels = c(
      total_body = "Total body length",
      total_leg = "Total leg length"
    )
  ) +
  labs(
    x = "Trait",
    y = "Trait length (mm)",
    fill = "Sex"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

ggsave("figure_2.png", width = 8, height = 8, dpi = 600)

summary_mating <- morph_data %>%
    group_by(sex, mated) %>%
    summarise(n = n())
  
summary_mating
# =============================
# selection analysis
# =============================

# male selection analysis
males <- morph_data %>%
  dplyr::filter(sex=="m")

# male.gam<-gam(mated~s(total_body)+s(total_leg),family=binomial(link="logit"), method="GCV.Cp", data=males)
# male.grad<-gam.gradients(mod = male.gam,phenotype=c("total_body","total_leg"), se.method="boot.para", n.boot=1000, standardize=TRUE)
# male.estimates<-data.frame(male.grad$ests)

## need to double the (quadratic) gamma gradient and SEs
# male.estimates[3,1]<-male.estimates[3,1]*2
# male.estimates[4,1]<-male.estimates[4,1]*2
# male.estimates[3,2]<-male.estimates[3,2]*2
# male.estimates[4,2]<-male.estimates[4,2]*2
# saveRDS(male.estimates, file = "data/processed/male.estimates.rds")

male.estimates <- readRDS("data/processed/male.estimates.rds")

# male selection differentials
# m.body <- moments.differentials(z=males$total_body, W=males$mated, n.boot=10000, standardized=TRUE)
# saveRDS(m.body, file = "data/processed/m.body.rds")
m.body <- readRDS("data/processed/m.body.rds")

# m.leg <- moments.differentials(z=males$total_leg, W=males$mated, n.boot=10000, standardized=TRUE)
# saveRDS(m.leg, file = "data/processed/m.leg.rds")
m.leg <- readRDS("data/processed/m.leg.rds")

# png("figure_4.png",
#     width = 2000, height = 2000, res = 300)
# 
# par(xpd = NA,
#     bg = "transparent",
#     oma = c(2, 2, 0, 0),
#     mgp = c(2, 1, 0))
# 
# vis.gam(male.gam,
#         view = c("total_body", "total_leg"),
#         color = "bw",
#         plot.type = "persp",
#         type = "response",
#         theta = 155,
#         xlab = "\nBody length (mm)",
#         ylab = "\nLeg length (mm)",
#         zlab = "\nFitness")
# 
# dev.off()

# female selection analysis
females <- morph_data %>%
  dplyr::filter(sex=="f")

# female.gam<-gam(mated~s(total_body)+s(total_leg),family=binomial(link="logit"),method="GCV.Cp", data=females)
# female.grad<-gam.gradients(mod = female.gam, phenotype=c("total_body","total_leg"),se.method="boot.para",n.boot=1000,standardize=TRUE)
# female.estimates<-data.frame(female.grad$ests)

# ## need to double the (quadratic) gamma gradient and SEs
# female.estimates[3,1]<-female.estimates[3,1]*2
# female.estimates[4,1]<-female.estimates[4,1]*2
# female.estimates[3,2]<-female.estimates[3,2]*2
# female.estimates[4,2]<-female.estimates[4,2]*2
# saveRDS(female.estimates, file = "data/processed/female.estimates.rds")

female.estimates <- readRDS("data/processed/female.estimates.rds")

# female selection differentials
# f.body <- moments.differentials(z=females$total_body,W=females$mated, n.boot=10000, standardized=TRUE)
# saveRDS(f.body, file = "data/processed/f.body.rds")
f.body <- readRDS("data/processed/f.body.rds")

# f.leg <- moments.differentials(z=females$total_leg,W=females$mated, n.boot=10000, standardized=TRUE)
# saveRDS(f.leg, file = "data/processed/f.leg.rds")
f.leg <- readRDS("data/processed/f.leg.rds")

# =========================
# Gradient plots
# =========================
# m.body.fitness <- fitness.landscape(mod=male.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
# saveRDS(m.body.fitness, file = "data/processed/m.body.fitness.rds")
m.body.fitness <- readRDS("data/processed/m.body.fitness.rds")

# m.leg.fitness <- fitness.landscape(mod=male.gam, phenotype="total_leg", covariates=c("total_body"), PI.method="boot.para")
# saveRDS(m.leg.fitness, file = "data/processed/m.leg.fitness.rds")
m.leg.fitness <- readRDS("data/processed/m.leg.fitness.rds")

# f.body.fitness <- fitness.landscape(mod=female.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
# saveRDS(f.body.fitness, file = "data/processed/f.body.fitness.rds")
f.body.fitness <- readRDS("data/processed/f.body.fitness.rds")

# f.leg.fitness <- fitness.landscape(mod=female.gam,phenotype="total_leg",covariates=c("total_body"),PI.method="boot.para")
# saveRDS(f.leg.fitness, file = "data/processed/f.leg.fitness.rds")
f.leg.fitness <- readRDS("data/processed/f.leg.fitness.rds")

# png("figure_1.png", width = 2400, height = 2400, res = 400, type = "cairo")
# 
# par(mfrow = c(2,2),
#     mar = c(5,5,3,1),
#     oma = c(0,0,2,0))
# 
# ### Helper function to avoid repetition
# plot_fitness <- function(obj, xlab, panel_label) {
#   
#   x <- obj$points[,1]
#   y <- obj$Wbar
#   lower <- obj$WbarPI[1,]
#   upper <- obj$WbarPI[2,]
#   
#   # Empty plot (sets axes)
#   plot(x, y, type="n",
#        ylim=c(0,0.8),
#        xlab=xlab,
#        ylab="Mating success",
#        cex.lab=1.6,
#        cex.axis=1.2,
#        bty = "l")
#   
#   # Shaded confidence interval
#   polygon(c(x, rev(x)),
#           c(upper, rev(lower)),
#           col=rgb(0,0,0,0.15),
#           border=NA)
#   
#   # Main line
#   lines(x, y, lwd=2.5)
#   
#   # CI lines
#   lines(x, lower, lty=2, lwd=1.2)
#   lines(x, upper, lty=2, lwd=1.2)
#   
#   # Panel label (top-left inside plot)
#   usr <- par("usr")
#   
#   text(x = usr[1] - 0.34 * diff(usr[1:2]),  # move left of y-axis label
#        y = usr[4] + 0.19 * diff(usr[3:4]),  # move ABOVE plot (and label)
#        labels = panel_label,
#        xpd = NA,
#        adj = c(0, 1),
#        cex = 1.4,
#        font = 2)
# }
# 
# ### Panels
# plot_fitness(m.body.fitness, "Body length (mm)", "(a)")
# plot_fitness(m.leg.fitness,  "Hind leg length (mm)", "(b)")
# plot_fitness(f.body.fitness, "Body length (mm)", "(c)")
# plot_fitness(f.leg.fitness,  "Hind leg length (mm)", "(d)")
# 
# dev.off()

# ===========================
# selection table
# ===========================

format_diff <- function(x, trait){
  x %>%
    rownames_to_column("component") %>%
    dplyr::select(component, Coefficient, SE, `P-value`) %>%
    mutate(trait = trait) %>%
    pivot_wider(
      names_from = component,
      values_from = c(Coefficient, SE, `P-value`),
      names_glue = "{.value}_{component}"
    ) %>%
    rename_with(~ gsub(" 1$", "", .x)) %>%
    rename_with(~ gsub("P-value", "pvalue", .x))
}

format_grad <- function(x){
  x %>%
    rownames_to_column("term") %>%
    mutate(
      type = case_when(
        str_detect(term, "^B-") ~ "beta",
        str_detect(term, "^G-") ~ "gamma"
      ),
      trait = str_remove(term, "^[BG]-"),
      trait = case_when(
        trait == "total_body" ~ "body length",
        trait == "total_leg"  ~ "leg length",
        TRUE ~ trait
      )
    ) %>%
    dplyr::select(trait, type, estimates, SE, P.value) %>%
    pivot_wider(
      names_from = type,
      values_from = c(estimates, SE, P.value),
      names_glue = "{type}_{.value}"
    ) %>%
    rename(
      beta     = beta_estimates,
      beta_se  = beta_SE,
      beta_p   = beta_P.value,
      gamma    = gamma_estimates,
      gamma_se = gamma_SE,
      gamma_p  = gamma_P.value
    )
}

male_table <- left_join(
  bind_rows(
    format_diff(m.body, "body length"),
    format_diff(m.leg,  "leg length")
  ),
  format_grad(male.estimates),
  by = "trait"
) %>%
  mutate(sex = "Male")

female_table <- left_join(
  bind_rows(
    format_diff(f.body, "body length"),
    format_diff(f.leg,  "leg length")
  ),
  format_grad(female.estimates),
  by = "trait"
) %>%
  mutate(sex = "Female")

final_table <- bind_rows(female_table, male_table) %>%
  mutate(
    S       = sprintf("%.3f ± %.3f", Coefficient_S, SE_S),
    S_p     = sprintf("%.3f", pvalue_S),
    C       = sprintf("%.3f ± %.3f", Coefficient_C, SE_C),
    C_p     = sprintf("%.3f", pvalue_C),
    beta    = sprintf("%.3f ± %.3f", beta, beta_se),
    beta_p  = sprintf("%.3f", beta_p),
    gamma   = sprintf("%.3f ± %.3f", gamma, gamma_se),
    gamma_p = sprintf("%.3f", gamma_p)
  ) %>%
  dplyr::select(sex, trait, S, S_p, C, C_p, beta, beta_p, gamma, gamma_p)

female_df <- final_table %>% filter(sex == "Female") %>% mutate(sex = "")
male_df   <- final_table %>% filter(sex == "Male")   %>% mutate(sex = "")

# helper: empty row
empty_row <- function(df) df[1, ] %>% mutate(across(everything(), ~ ""))

# section headers
female_header <- empty_row(final_table)
female_header$sex <- "(a) Females"

male_header <- empty_row(final_table)
male_header$sex <- "(b) Males"

# spacer
spacer <- empty_row(final_table)
spacer[] <- ""

table_combined <- bind_rows(
  female_header,
  female_df,
  spacer,
  male_header,
  male_df
)

table_combined <- table_combined %>%
  mutate(
    P_S_num = suppressWarnings(as.numeric(S_p)),
    P_C_num = suppressWarnings(as.numeric(C_p)),
    P_β_num = suppressWarnings(as.numeric(beta_p)),
    P_γ_num = suppressWarnings(as.numeric(gamma_p))
  )

table_combined <- table_combined %>%
  mutate(
    P_S_fmt = ifelse(is.na(P_S_num), "", ifelse(P_S_num < 0.001, "<0.001", sprintf("%.3f", P_S_num))),
    P_C_fmt = ifelse(is.na(P_C_num), "", ifelse(P_C_num < 0.001, "<0.001", sprintf("%.3f", P_C_num))),
    P_β_fmt = ifelse(is.na(P_β_num), "", ifelse(P_β_num < 0.001, "<0.001", sprintf("%.3f", P_β_num))),
    P_γ_fmt = ifelse(is.na(P_γ_num), "", ifelse(P_γ_num < 0.001, "<0.001", sprintf("%.3f", P_γ_num)))
  )

table_combined <- table_combined %>%
  rename(
    Sex   = sex,
    Trait = trait,
    `S ± SE` = S,
    `P_S`    = P_S_fmt,
    `C ± SE` = C,
    `P_C`    = P_C_fmt,
    `β ± SE` = beta,
    `P_β`    = P_β_fmt,
    `γ ± SE` = gamma,
    `P_γ`    = P_γ_fmt
  )

table_combined <- table_combined %>%
  mutate(spacer = "") %>%
  relocate(spacer, .after = `P_C`)

table_combined <- table_combined %>%
  dplyr::select(
    Sex, Trait,
    `S ± SE`, P_S,
    `C ± SE`, P_C,
    spacer,
    `β ± SE`, P_β,
    `γ ± SE`, P_γ,
    P_S_num, P_C_num, P_β_num, P_γ_num
  )

table_combined_display <- table_combined %>%
  dplyr::select(
    Sex, Trait,
    `S ± SE`, P_S,
    `C ± SE`, P_C,
    spacer,
    `β ± SE`, P_β,
    `γ ± SE`, P_γ
  )

ft <- flextable(table_combined_display)

cols <- ft$col_keys
spacer_col <- which(cols == "spacer")

ft <- add_header_row(
  ft,
  values = c("", "", "Differentials", "", "Gradients"),
  colwidths = c(1, 1, 4, 1, 4)
)

ft <- set_header_labels(
  ft,
  `P_S` = "P",
  `P_C` = "P",
  `P_β` = "P",
  `P_γ` = "P",
  spacer = "",
  Sex = ""
)

ft <- border_remove(ft)

left_cols  <- seq_len(spacer_col - 1)
right_cols <- (spacer_col + 1):length(cols)

# Row 1 (split line)
# define columns EXCLUDING spacer
diff_cols <- which(cols %in% c("S ± SE", "P_S", "C ± SE", "P_C"))
grad_cols <- which(cols %in% c("β ± SE", "P_β", "γ ± SE", "P_γ"))

# ensure spacer is excluded (critical)
diff_cols <- setdiff(diff_cols, spacer_col)
grad_cols <- setdiff(grad_cols, spacer_col)

# draw split lines
ft <- border(ft, i = 1, j = diff_cols,
             border.bottom = fp_border(width = 1), part = "header")

ft <- border(ft, i = 1, j = grad_cols,
             border.bottom = fp_border(width = 1), part = "header")

# Row 2 (FULL continuous line)
ft <- border(ft, i = 2, j = seq_along(cols),
             border.bottom = fp_border(width = 1.5), part = "header")

ft <- border(
  ft,
  j = spacer_col,
  border.left  = fp_border(width = 0),
  border.right = fp_border(width = 0),
  part = "all"
)

section_rows <- which(table_combined$Sex %in% c("(a) Females", "(b) Males"))

for (i in section_rows) {
  ft <- merge_at(ft, i = i, j = seq_along(cols))
}

# general alignment
ft <- ft %>%
  align(part = "header", align = "center") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  align(j = setdiff(seq_along(cols), spacer_col),
        align = "center", part = "all")

# left-align section headers
ft <- align(
  ft,
  i = section_rows,
  j = seq_along(cols),
  align = "left",
  part = "body"
)

trait_rows <- setdiff(seq_len(nrow(table_combined)), section_rows)

ft <- padding(
  ft,
  i = trait_rows,
  j = "Trait",
  padding.left = 10
)

ft <- ft %>%
  bold(part = "header") %>%
  bold(i = section_rows, bold = TRUE) %>%
  autofit()

ft <- border(ft, i = nrow(table_combined), j = seq_along(cols),
             border.bottom = fp_border(width = 1.5), part = "body")

ft <- border(
  ft,
  i = 1,
  j = seq_along(cols),
  border.top = fp_border(width = 1.5),
  part = "header"
)

# shrink spacer column AFTER autofit
ft <- width(ft, j = spacer_col, width = 0.05)

alpha <- 0.05

ft <- bold(ft,
           i = which(table_combined$P_S_num < alpha),
           j = "S ± SE",
           part = "body"
)

ft <- bold(ft,
           i = which(table_combined$P_C_num < alpha),
           j = "C ± SE",
           part = "body"
)

ft <- bold(ft,
           i = which(table_combined$P_β_num < alpha),
           j = "β ± SE",
           part = "body"
)

ft <- bold(ft,
           i = which(table_combined$P_γ_num < alpha),
           j = "γ ± SE",
           part = "body"
)

ft <- bold(ft, i = ~ P_S < alpha, j = "P_S", part = "body")
ft <- bold(ft, i = ~ P_C < alpha, j = "P_C", part = "body")
ft <- bold(ft, i = ~ P_β < alpha, j = "P_β", part = "body")
ft <- bold(ft, i = ~ P_γ < alpha, j = "P_γ", part = "body")

ft

# =========================
# female fecundity
# =========================

females_clean <- females |>
  filter(!is.na(eggs), total_body > 0)

female.fecundity <- glm.nb(eggs ~ log(total_body), data=females_clean)
summary(female.fecundity) # strong fecundity selection

# newdat <- data.frame(
#   total_body = seq(min(females$total_body),
#                  max(females$total_body),
#                  length.out = 100)
# )
# 
# pred <- predict(female.fecundity, newdata = newdat,
#                 type = "link", se.fit = TRUE)
# 
# newdat$fit <- exp(pred$fit)
# newdat$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
# newdat$upr <- exp(pred$fit + 1.96 * pred$se.fit)
# 
# sum(is.na(females$total_body))
# sum(is.na(females$eggs))
# 
# sum(females$total_body <= 0, na.rm = TRUE)
# 
# female.fecundity.plot <- ggplot(females_clean, aes(x = total_body, y = eggs)) +
#       geom_point(alpha = 0.5) +
#       geom_line(data = newdat, aes(x = total_body, y = fit),
#             color = "blue", linewidth = 1) +
#       geom_ribbon(data = newdat,
#               aes(x = total_body, ymin = lwr, ymax = upr),
#               inherit.aes = FALSE,
#               alpha = 0.2) +
#       scale_x_log10() +
#       labs(x = "Body size (mm, log scale)",
#                y = "Fecundity (number of eggs)") +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = c(0.98, 0.98),
#     legend.justification = c(1, 1),
#     legend.background = element_blank(),
#     legend.key = element_blank(),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(color = "black")
#   )
# 
# ggsave("figure_3.png", width = 8, height = 8, dpi = 600)


# =========================
# assortative mating
# =========================

pairs <- morph_data %>%
  filter(pair!="NA")%>%
  pivot_wider(
    id_cols = pair,
    names_from = sex,
    values_from = c(total_body, total_leg, mated)
  )

assort <- lm(scale(total_body_f) ~ scale(total_body_m), data = pairs)
summary(assort)

obs_slope <- coef(assort)[2]

null_slopes <- replicate(10000, {
  perm_females <- sample(pairs$total_body_f)
  coef(lm(perm_females ~ pairs$total_body_m))[2]
})

p_val <- mean(abs(null_slopes) >= abs(obs_slope)) # no assortative mating

## ---- end

# =========================
# male mate choice
# =========================

# male_mate_choice <- read.csv(file="data/raw/male_mate_choice.csv", header=TRUE, sep=",", dec=".") %>%
#   as.data.frame()
# 
# choice <- male_mate_choice %>%
#   left_join(
#     females %>% dplyr::select(id, total_body, total_leg),
#     by = c("chosen_female_id" = "id")
#   ) %>%
#   rename(
#     total_body_chosen = total_body,
#     leg_chosen        = total_leg
#   ) %>%
#   left_join(
#     females %>% dplyr::select(id, total_body, total_leg),
#     by = c("unchosen_fem_id" = "id")
#   ) %>%
#   rename(
#     total_body_unchosen = total_body,
#     leg_unchosen        = total_leg
#   )
# 
# choice$chose_larger <- ifelse(choice$total_body_chosen > choice$total_body_unchosen, 1, 0)
# choice$size_diff <- choice$total_body_chosen - choice$total_body_unchosen
# 
# table(choice$chose_larger)
# 
# binom.test(sum(choice$chose_larger), nrow(choice), p = 0.5)

# =========================
# male-male competition
# =========================

male_male <- read.csv(file="data/raw/male_male_combat.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

# dyadic
combat2 <- male_male %>%
  left_join(
    males %>% dplyr::select(id, total_body, total_leg),
    by = c("winner_id" = "id")
  ) %>%
  rename(
    total_body_win = total_body,
    leg_win        = total_leg
  ) %>%
  left_join(
    males %>% dplyr::select(id, total_body, total_leg),
    by = c("loser_id" = "id")
  ) %>%
  rename(
    total_body_lose = total_body,
    leg_lose        = total_leg
  )

# combat2$won_larger <- ifelse(combat2$total_body_win > combat2$total_body_lose, 1, 0)
# combat2$size_diff <- combat2$total_body_win - combat2$total_body_lose
# 
# binom.test(sum(combat2$won_larger), nrow(combat2), p = 0.5)

combat2 <- combat2 %>%
  mutate(
    body_diff = total_body_win - total_body_lose,
    leg_diff  = leg_win - leg_lose
  ) # positive = winner is bigger

df <- combat2

set.seed(123)

flip <- rbinom(nrow(df), 1, 0.5)

df$id_i <- ifelse(flip == 1, df$winner_id, df$loser_id)
df$id_j <- ifelse(flip == 1, df$loser_id, df$winner_id)

df$leg_i <- ifelse(flip == 1, df$leg_win, df$leg_lose)
df$leg_j <- ifelse(flip == 1, df$leg_lose, df$leg_win)

df$body_i <- ifelse(flip == 1, df$total_body_win, df$total_body_lose)
df$body_j <- ifelse(flip == 1, df$total_body_lose, df$total_body_win)

df$win <- ifelse(df$id_i == df$winner_id, 1, 0)

df$diff_leg  <- df$leg_i - df$leg_j
df$diff_body <- df$body_i - df$body_j

model.dyad <- glm(win ~ scale(diff_leg) + scale(diff_body),
             family = binomial, data = df)

summary(model.dyad)

# ggplot(df, aes(diff_body, win)) +
#   geom_jitter(height = 0.05) +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))

## reaction norm plot

df_long <- df %>%
 pivot_longer(
    cols = starts_with("total_body"),
    names_to = "status",
    values_to = "total_body"
  ) %>%
  mutate(
    status = ifelse(status == "total_body_win", "winner", "loser")
  ) %>%
  print(n=100)

ggplot(df_long, aes(x = status, y = total_body, group = trial)) +
  geom_line(alpha = 0.5) +
  geom_point() +
  theme_classic()

## selection analysis
winners <- data.frame(
  id   = combat2$winner_id,
  leg  = combat2$leg_win,
  body = combat2$total_body_win,
  win  = 1
)

losers <- data.frame(
  id   = combat2$loser_id,
  leg  = combat2$leg_lose,
  body = combat2$total_body_lose,
  win  = 0
)

indiv <- rbind(winners, losers)

indiv$leg_z  <- scale(indiv$leg)
indiv$body_z <- scale(indiv$body)

model_sel <- glm(win ~ leg_z + body_z, family = binomial, data=indiv)
summary(model_sel)







