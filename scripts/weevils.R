# If Git is ahead by X commits, use: git reset --soft HEAD~X

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
library(broman)
library(scales)
# install.packages("~/Documents/Action/gsg_2.0.tar", repos = NULL, type = "source")

# Create the output directory if needed
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

#### read in data file ####
# Load raw morphometric data and calculate derived traits
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
# Ensure sex is treated as a categorical factor for modeling
morph_data$sex <- factor(morph_data$sex)

# Sex-specific datasets used in the selection analyses
males <- morph_data %>%
  filter(sex == "m") %>%
  drop_na(mated, total_body, total_leg)

females <- morph_data %>%
  filter(sex == "f") %>%
  drop_na(mated, total_body, total_leg)

# Confirm the sample sizes used in the pairing-status analyses
print(males %>% count(mated))
print(females %>% count(mated))

stopifnot(
  nrow(males) == 100,
  nrow(females) == 77,
  sum(males$mated == 1) == 25,
  sum(females$mated == 1) == 25
)

str(morph_data)
# Inspect unusually short male legs, if present
morph_data %>%
  filter(sex == "m", total_leg < 6)

# Calculate summary statistics by sex
summary_stats <- morph_data %>%
  group_by(sex) %>%
  summarise(
    n_body = sum(!is.na(total_body)),
    mean_total_body = mean(total_body, na.rm = TRUE),
    se_total_body = sd(total_body, na.rm = TRUE) / sqrt(n_body),
    
    n_leg = sum(!is.na(total_leg)),
    mean_total_leg = mean(total_leg, na.rm = TRUE),
    se_total_leg = sd(total_leg, na.rm = TRUE) / sqrt(n_leg),
    .groups = "drop"
  )

summary_stats

# total body length
# Fit linear models to test for sex differences in body and leg length
model_body <- lm(total_body ~ sex, data = morph_data)
summary(model_body)

# total leg length
model_leg <- lm(total_leg ~ sex, data = morph_data)
summary(model_leg)

# Test for multivariate sexual dimorphism using MANOVA
manova_model <- manova(cbind(total_body, total_leg) ~ sex, data = morph_data)
summary(manova_model)

# Test whether the static allometric relationship between hind leg length and
# body length differs between the sexes. Log transformation makes the slope an
# allometric exponent and the sex-by-body interaction tests slope heterogeneity.
allometry_model <- lm(
  log(total_leg) ~ log(total_body) * sex,
  data = morph_data
)
summary(allometry_model)

plot_data <- morph_data %>%
  dplyr::select(sex, total_body, total_leg) %>%
  pivot_longer(
    cols = c(total_body, total_leg),
    names_to = "trait",
    values_to = "value"
  )

# figure_5 <- ggplot(plot_data, aes(x = trait, y = value, fill = sex)) +
#   geom_boxplot(
#     position = position_dodge(width = 0.7),
#     width = 0.6,
#     outlier.shape = NA,
#     color = "black"
#   ) +
#   geom_point(
#     position = position_jitterdodge(
#       jitter.width = 0.15,
#       dodge.width = 0.7
#     ),
#     size = 1.6,
#     alpha = 0.7,
#     color = "black"
#   ) +
#   scale_fill_grey(
#     start = 0.8,
#     end = 0.4,
#     labels = c("f" = "Females", "m" = "Males")
#   ) +
#   scale_x_discrete(
#     labels = c(
#       total_body = "Total body length",
#       total_leg = "Total leg length"
#     )
#   ) +
#   labs(
#     x = "Trait",
#     y = "Trait length (mm)",
#     fill = "Sex"
#   ) +
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
# ggsave("figure_5.png", width = 8, height = 8, dpi = 600)

# Summarize mating status by sex
summary_mating <- morph_data %>%
    group_by(sex, mated) %>%
    summarise(n = n())
  
summary_mating
# =============================
# selection analysis
# =============================

# Selection differentials: 10,000 nonparametric bootstrap replicates
set.seed(2025)

# m.body <- moments.differentials(
#   z = males["total_body"],
#   W = males$mated,
#   n.boot = 10000,
#   standardized = TRUE
# )
# 
# m.leg <- moments.differentials(
#   z = males["total_leg"],
#   W = males$mated,
#   n.boot = 10000,
#   standardized = TRUE
# )
# 
# f.body <- moments.differentials(
#   z = females["total_body"],
#   W = females$mated,
#   n.boot = 10000,
#   standardized = TRUE
# )
# 
# f.leg <- moments.differentials(
#   z = females["total_leg"],
#   W = females$mated,
#   n.boot = 10000,
#   standardized = TRUE
# )

# saveRDS(m.body, "data/processed/m.body.rds")
# saveRDS(m.leg, "data/processed/m.leg.rds")
# saveRDS(f.body, "data/processed/f.body.rds")
# saveRDS(f.leg, "data/processed/f.leg.rds")

m.body <- readRDS("data/processed/m.body.rds")
m.leg <- readRDS("data/processed/m.leg.rds")
f.body <- readRDS("data/processed/f.body.rds")
f.leg <- readRDS("data/processed/f.leg.rds")

# Explicit GAM specifications
male.gam <- gam(
  mated ~
    s(total_body, bs = "tp", k = 10) +
    s(total_leg,  bs = "tp", k = 10),
  family = binomial(link = "logit"),
  method = "GCV.Cp",
  data = males
)

female.gam <- gam(
  mated ~
    s(total_body, bs = "tp", k = 10) +
    s(total_leg,  bs = "tp", k = 10),
  family = binomial(link = "logit"),
  method = "GCV.Cp",
  data = females
)

# Selection gradients: 10,000 parametric-bootstrap replicates
set.seed(2026)

# male.grad <- gam.gradients(
#   mod = male.gam,
#   phenotype = c("total_body", "total_leg"),
#   se.method = "boot.para",
#   n.boot = 10000,
#   standardized = TRUE,
#   refit.smooth = FALSE
# )

# female.grad <- gam.gradients(
#   mod = female.gam,
#   phenotype = c("total_body", "total_leg"),
#   se.method = "boot.para",
#   n.boot = 10000,
#   standardized = TRUE,
#   refit.smooth = FALSE
# )

# male.estimates <- as.data.frame(male.grad$ests)
# female.estimates <- as.data.frame(female.grad$ests)
# 
# saveRDS(male.estimates, "data/processed/male.estimates.rds")
# saveRDS(female.estimates, "data/processed/female.estimates.rds")

male.estimates <- readRDS("data/processed/male.estimates.rds")
female.estimates <- readRDS("data/processed/female.estimates.rds")

# Figure 3: explicit 95% parametric-bootstrap intervals
set.seed(2026)

# m.body.fitness <- fitness.landscape(
#   mod = male.gam,
#   phenotype = "total_body",
#   covariates = "total_leg",
#   plt.density = 25,
#   PI.method = "boot.para",
#   PI.interval = c(0.025, 0.975),
#   n.boot = 1000,
#   refit.smooth = FALSE
# )
# 
# m.leg.fitness <- fitness.landscape(
#   mod = male.gam,
#   phenotype = "total_leg",
#   covariates = "total_body",
#   plt.density = 25,
#   PI.method = "boot.para",
#   PI.interval = c(0.025, 0.975),
#   n.boot = 1000,
#   refit.smooth = FALSE
# )
# 
# f.body.fitness <- fitness.landscape(
#   mod = female.gam,
#   phenotype = "total_body",
#   covariates = "total_leg",
#   plt.density = 25,
#   PI.method = "boot.para",
#   PI.interval = c(0.025, 0.975),
#   n.boot = 1000,
#   refit.smooth = FALSE
# )
# 
# f.leg.fitness <- fitness.landscape(
#   mod = female.gam,
#   phenotype = "total_leg",
#   covariates = "total_body",
#   plt.density = 25,
#   PI.method = "boot.para",
#   PI.interval = c(0.025, 0.975),
#   n.boot = 1000,
#   refit.smooth = FALSE
# )
# 
# saveRDS(m.body.fitness, "data/processed/m.body.fitness.rds")
# saveRDS(m.leg.fitness,  "data/processed/m.leg.fitness.rds")
# saveRDS(f.body.fitness, "data/processed/f.body.fitness.rds")
# saveRDS(f.leg.fitness,  "data/processed/f.leg.fitness.rds")
# 
# #----------------------------------------------------------
# # 1. Extract the observations used in the model
# #----------------------------------------------------------

# male_data <- model.frame(male.gam)
# 
# male_data <- male_data[
#   complete.cases(
#     male_data$total_body,
#     male_data$total_leg
#   ),
# ]
# 
# #----------------------------------------------------------
# # 2. Construct a prediction grid
# #----------------------------------------------------------
# 
# x_seq <- seq(
#   min(male_data$total_body),
#   max(male_data$total_body),
#   length.out = 300
# )
# 
# y_seq <- seq(
#   min(male_data$total_leg),
#   max(male_data$total_leg),
#   length.out = 300
# )
# 
# pred_grid <- expand.grid(
#   total_body = x_seq,
#   total_leg  = y_seq
# )
# 
# pred_grid$probability <- predict(
#   male.gam,
#   newdata = pred_grid,
#   type = "response"
# )
# 
# #----------------------------------------------------------
# # 3. Identify the convex hull of the observations
# #----------------------------------------------------------
# 
# hull_id <- chull(
#   male_data$total_body,
#   male_data$total_leg
# )
# 
# hull_x <- male_data$total_body[hull_id]
# hull_y <- male_data$total_leg[hull_id]
# 
# # Close the polygon
# hull_x <- c(hull_x, hull_x[1])
# hull_y <- c(hull_y, hull_y[1])
# 
# # Determine which prediction cells fall within the hull
# inside_hull <- sp::point.in.polygon(
#   point.x = pred_grid$total_body,
#   point.y = pred_grid$total_leg,
#   pol.x   = hull_x,
#   pol.y   = hull_y
# ) > 0
# 
# # Mask unsupported trait combinations
# pred_grid$probability[!inside_hull] <- NA_real_
# 
# # Convert predictions to the matrix required by image()
# z_matrix <- matrix(
#   pred_grid$probability,
#   nrow = length(x_seq),
#   ncol = length(y_seq)
# )
# 
# #----------------------------------------------------------
# # 4. Draw the figure
# #----------------------------------------------------------
# 
# png(
#   "figure_4.png",
#   width = 2200,
#   height = 1900,
#   res = 300,
#   bg = "white"
# )
# 
# par(
#   mar = c(5, 5, 1, 1),
#   mgp = c(3, 1, 0),
#   xaxs = "i",
#   yaxs = "i"
# )
# 
# image(
#   x = x_seq,
#   y = y_seq,
#   z = z_matrix,
#   col = hcl.colors(100, palette = "Viridis"),
#   zlim = c(0, 1),
#   useRaster = TRUE,
#   xlab = "Body length (mm)",
#   ylab = "Leg length (mm)",
#   main = ""
# )
# 
# contour(
#   x = x_seq,
#   y = y_seq,
#   z = z_matrix,
#   levels = seq(0.1, 0.9, by = 0.1),
#   add = TRUE,
#   drawlabels = TRUE,
#   labcex = 0.75,
#   lwd = 0.8,
#   col = "black"
# )
# 
# points(
#   male_data$total_body,
#   male_data$total_leg,
#   pch = 16,
#   cex = 0.4,
#   col = adjustcolor("black", alpha.f = 0.4)
# )
# 
# box()
# 
# dev.off()

# Regenerate Figure 4 with high-contrast symbols so observations remain
# visible across the full colour surface. Pairing status is also encoded by
# symbol and fill rather than colour alone.
make_figure_4 <- function(
    model = male.gam,
    output_file = "figure_4.png"
) {
  male_data <- model.frame(model)
  male_data <- male_data[
    complete.cases(male_data$total_body, male_data$total_leg, male_data$mated),
  ]

  x_seq <- seq(
    min(male_data$total_body),
    max(male_data$total_body),
    length.out = 300
  )
  y_seq <- seq(
    min(male_data$total_leg),
    max(male_data$total_leg),
    length.out = 300
  )

  pred_grid <- expand.grid(
    total_body = x_seq,
    total_leg = y_seq
  )
  pred_grid$probability <- predict(
    model,
    newdata = pred_grid,
    type = "response"
  )

  hull_id <- chull(male_data$total_body, male_data$total_leg)
  hull_x <- male_data$total_body[hull_id]
  hull_y <- male_data$total_leg[hull_id]

  # Ray-casting point-in-polygon test avoids adding a plotting dependency.
  inside_hull <- rep(FALSE, nrow(pred_grid))
  previous_vertex <- length(hull_x)
  for (current_vertex in seq_along(hull_x)) {
    edge_crosses <-
      (hull_y[current_vertex] > pred_grid$total_leg) !=
      (hull_y[previous_vertex] > pred_grid$total_leg)
    crossing_x <-
      (hull_x[previous_vertex] - hull_x[current_vertex]) *
      (pred_grid$total_leg - hull_y[current_vertex]) /
      (hull_y[previous_vertex] - hull_y[current_vertex]) +
      hull_x[current_vertex]
    inside_hull <- xor(
      inside_hull,
      edge_crosses & pred_grid$total_body < crossing_x
    )
    previous_vertex <- current_vertex
  }
  pred_grid$probability[!inside_hull] <- NA_real_

  z_matrix <- matrix(
    pred_grid$probability,
    nrow = length(x_seq),
    ncol = length(y_seq)
  )

  png(
    output_file,
    width = 2200,
    height = 1900,
    res = 300,
    bg = "white"
  )
  on.exit(dev.off(), add = TRUE)

  par(
    mar = c(5, 5, 1, 1),
    mgp = c(3, 1, 0),
    xaxs = "i",
    yaxs = "i"
  )

  image(
    x = x_seq,
    y = y_seq,
    z = z_matrix,
    col = hcl.colors(100, palette = "Viridis"),
    zlim = c(0, 1),
    useRaster = TRUE,
    xlab = "Body length (mm)",
    ylab = "Leg length (mm)"
  )
  contour(
    x = x_seq,
    y = y_seq,
    z = z_matrix,
    levels = seq(0.1, 0.9, by = 0.1),
    add = TRUE,
    drawlabels = TRUE,
    labcex = 0.75,
    lwd = 0.8,
    col = "black"
  )

  unpaired <- male_data$mated == 0
  points(
    male_data$total_body[unpaired],
    male_data$total_leg[unpaired],
    pch = 21,
    cex = 0.70,
    lwd = 0.9,
    col = "black",
    bg = "white"
  )
  points(
    male_data$total_body[!unpaired],
    male_data$total_leg[!unpaired],
    pch = 24,
    cex = 0.85,
    lwd = 1.0,
    col = "black",
    bg = "#FDE725"
  )
  legend(
    "topleft",
    inset = 0.015,
    legend = c("Unpaired", "Paired"),
    pch = c(21, 24),
    pt.bg = c("white", "#FDE725"),
    col = "black",
    pt.cex = c(0.9, 1.0),
    bty = "n"
  )
  box()
}

make_figure_4()

# =========================
# Gradient plots
# =========================
# m.body.fitness <- fitness.landscape(mod=male.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
# saveRDS(m.body.fitness, file = "data/processed/m.body.fitness.rds")
# Load fitness landscape data for plotting gradients
m.body.fitness <- readRDS("data/processed/m.body.fitness.rds")

# m.leg.fitness <- fitness.landscape(mod=male.gam, phenotype="total_leg", covariates=c("total_body"), PI.method="boot.para")
# saveRDS(m.leg.fitness, file = "data/processed/m.leg.fitness.rds")
m.leg.fitness <- readRDS("data/processed/m.leg.fitness.rds")

# f.body.fitness <- fitness.landscape(mod=female.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
# saveRDS(f.body.fitness, file = "data/processed/f.body.fitness.rds")
# Load female fitness landscape data
f.body.fitness <- readRDS("data/processed/f.body.fitness.rds")

# f.leg.fitness <- fitness.landscape(mod=female.gam,phenotype="total_leg",covariates=c("total_body"),PI.method="boot.para")
# saveRDS(f.leg.fitness, file = "data/processed/f.leg.fitness.rds")
f.leg.fitness <- readRDS("data/processed/f.leg.fitness.rds")

png("figure_3.png", width = 2400, height = 2400, res = 400, type = "cairo")

par(mfrow = c(2,2),
    mar = c(5,5,3,1),
    oma = c(0,0,2,0))

### Helper function to avoid repetition
plot_fitness <- function(obj, xlab, panel_label) {

  x <- obj$points[,1]
  y <- obj$Wbar
  lower <- obj$WbarPI[1,]
  upper <- obj$WbarPI[2,]

  # Empty plot (sets axes)
  plot(x, y, type="n",
       ylim = c(0, 1),
       xlab=xlab,
       ylab="Predicted probability of \nbeing observed paired",
       cex.lab=1.0,
       cex.axis=1.0,
       bty = "l")

  # Shaded confidence interval
  polygon(c(x, rev(x)),
          c(upper, rev(lower)),
          col=rgb(0,0,0,0.15),
          border=NA)

  # Main line
  lines(x, y, lwd=2.5)

  # CI lines
  lines(x, lower, lty=2, lwd=1.2)
  lines(x, upper, lty=2, lwd=1.2)

  # Panel label (top-left inside plot)
  usr <- par("usr")

  text(x = usr[1] - 0.34 * diff(usr[1:2]),  # move left of y-axis label
       y = usr[4] + 0.19 * diff(usr[3:4]),  # move ABOVE plot (and label)
       labels = panel_label,
       xpd = NA,
       adj = c(0, 1),
       cex = 1.4,
       font = 2)
}

### Panels
plot_fitness(m.body.fitness, "Body length (mm)", "(a)")
plot_fitness(m.leg.fitness,  "Hind leg length (mm)", "(b)")
plot_fitness(f.body.fitness, "Body length (mm)", "(c)")
plot_fitness(f.leg.fitness,  "Hind leg length (mm)", "(d)")

dev.off()

# ===========================
# Selection table
# ===========================

#----------------------------------------------------------
# Helper function: format selection differentials
#----------------------------------------------------------

format_diff <- function(x, trait_name) {
  
  x %>%
    rownames_to_column("component") %>%
    dplyr::select(
      component,
      Coefficient,
      SE,
      `P-value`
    ) %>%
    mutate(trait = trait_name) %>%
    pivot_wider(
      names_from = component,
      values_from = c(Coefficient, SE, `P-value`),
      names_glue = "{.value}_{component}"
    ) %>%
    rename_with(~ gsub(" 1$", "", .x)) %>%
    rename_with(~ gsub("P-value", "pvalue", .x))
}

#----------------------------------------------------------
# Helper function: format selection gradients
#----------------------------------------------------------

format_grad <- function(x) {
  
  x %>%
    rownames_to_column("term") %>%
    mutate(
      type = case_when(
        str_detect(term, "^B-") ~ "beta",
        str_detect(term, "^G-") ~ "gamma",
        TRUE                    ~ NA_character_
      ),
      
      trait = str_remove(term, "^[BG]-"),
      
      trait = case_when(
        trait == "total_body" ~ "body length",
        trait == "total_leg" ~ "hind leg length",
        trait == "total_body-total_leg" ~
          "body length × hind leg length",
        TRUE ~ trait
      )
    ) %>%
    filter(!is.na(type)) %>%
    dplyr::select(
      trait,
      type,
      estimates,
      SE,
      P.value
    ) %>%
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

#----------------------------------------------------------
# Combine selection differentials and gradients
#----------------------------------------------------------

male_table <- full_join(
  bind_rows(
    format_diff(m.body, "body length"),
    format_diff(m.leg,  "hind leg length")
  ),
  format_grad(male.estimates),
  by = "trait"
) %>%
  mutate(
    sex = "Male",
    trait = factor(
      trait,
      levels = c(
        "body length",
        "hind leg length",
        "body length × hind leg length"
      )
    )
  ) %>%
  arrange(trait) %>%
  mutate(trait = as.character(trait))

female_table <- full_join(
  bind_rows(
    format_diff(f.body, "body length"),
    format_diff(f.leg,  "hind leg length")
  ),
  format_grad(female.estimates),
  by = "trait"
) %>%
  mutate(
    sex = "Female",
    trait = factor(
      trait,
      levels = c(
        "body length",
        "hind leg length",
        "body length × hind leg length"
      )
    )
  ) %>%
  arrange(trait) %>%
  mutate(trait = as.character(trait))

#----------------------------------------------------------
# Format estimates and standard errors
#----------------------------------------------------------

format_est_se <- function(estimate, se) {
  output <- sprintf("%.3f ± %.3f", estimate, se)
  output[is.na(estimate) | is.na(se)] <- ""
  output
}

format_raw_p <- function(x) {
  output <- sprintf("%.4f", x)
  output[is.na(x)] <- ""
  output
}

final_table <- bind_rows(
  female_table,
  male_table
) %>%
  mutate(
    S       = format_est_se(Coefficient_S, SE_S),
    S_p     = format_raw_p(pvalue_S),
    
    C       = format_est_se(Coefficient_C, SE_C),
    C_p     = format_raw_p(pvalue_C),
    
    beta    = format_est_se(beta, beta_se),
    beta_p  = format_raw_p(beta_p),
    
    gamma   = format_est_se(gamma, gamma_se),
    gamma_p = format_raw_p(gamma_p)
  ) %>%
  dplyr::select(
    sex,
    trait,
    S,
    S_p,
    C,
    C_p,
    beta,
    beta_p,
    gamma,
    gamma_p
  )

#----------------------------------------------------------
# Separate female and male rows
#----------------------------------------------------------

female_df <- final_table %>%
  filter(sex == "Female") %>%
  mutate(sex = "")

male_df <- final_table %>%
  filter(sex == "Male") %>%
  mutate(sex = "")

#----------------------------------------------------------
# Helper for creating blank rows
#----------------------------------------------------------

empty_row <- function(df) {
  
  df[1, , drop = FALSE] %>%
    mutate(across(everything(), ~ ""))
}

# Section-heading rows
female_header <- empty_row(final_table)
female_header$sex <- "(a) Females"

male_header <- empty_row(final_table)
male_header$sex <- "(b) Males"

# Blank spacer row between sections
spacer <- empty_row(final_table)

# Combine all table rows
table_combined <- bind_rows(
  female_header,
  female_df,
  spacer,
  male_header,
  male_df
)

#----------------------------------------------------------
# Retain numeric P-values for conditional formatting
#----------------------------------------------------------

table_combined <- table_combined %>%
  mutate(
    P_S_num = suppressWarnings(as.numeric(S_p)),
    P_C_num = suppressWarnings(as.numeric(C_p)),
    P_β_num = suppressWarnings(as.numeric(beta_p)),
    P_γ_num = suppressWarnings(as.numeric(gamma_p))
  )

#----------------------------------------------------------
# Format P-values for display
#----------------------------------------------------------

format_p <- function(x) {
  
  case_when(
    is.na(x)  ~ "",
    x < 0.001 ~ "<0.001",
    TRUE      ~ sprintf("%.3f", x)
  )
}

table_combined <- table_combined %>%
  mutate(
    P_S_fmt = format_p(P_S_num),
    P_C_fmt = format_p(P_C_num),
    P_β_fmt = format_p(P_β_num),
    P_γ_fmt = format_p(P_γ_num)
  )

#----------------------------------------------------------
# Rename displayed columns
#----------------------------------------------------------

table_combined <- table_combined %>%
  rename(
    Sex      = sex,
    Trait    = trait,
    `S ± SE` = S,
    P_S      = P_S_fmt,
    `C ± SE` = C,
    P_C      = P_C_fmt,
    `β ± SE` = beta,
    P_β      = P_β_fmt,
    `γ ± SE` = gamma,
    P_γ      = P_γ_fmt
  )

# Add narrow spacer column between the two table sections
table_combined <- table_combined %>%
  mutate(spacer = "") %>%
  relocate(spacer, .after = P_C)

# Arrange columns, retaining numeric P-values at the end
table_combined <- table_combined %>%
  dplyr::select(
    Sex,
    Trait,
    `S ± SE`,
    P_S,
    `C ± SE`,
    P_C,
    spacer,
    `β ± SE`,
    P_β,
    `γ ± SE`,
    P_γ,
    P_S_num,
    P_C_num,
    P_β_num,
    P_γ_num
  )

# Data used to construct the visible flextable
table_combined_display <- table_combined %>%
  dplyr::select(
    Sex,
    Trait,
    `S ± SE`,
    P_S,
    `C ± SE`,
    P_C,
    spacer,
    `β ± SE`,
    P_β,
    `γ ± SE`,
    P_γ
  )

#----------------------------------------------------------
# Create the flextable
#----------------------------------------------------------

ft <- flextable(table_combined_display)

cols <- ft$col_keys
spacer_col <- which(cols == "spacer")

#----------------------------------------------------------
# Add spanning header row
#----------------------------------------------------------

ft <- add_header_row(
  ft,
  values = c(
    "",
    "",
    "Differentials",
    "",
    "Gradients"
  ),
  colwidths = c(
    1,
    1,
    4,
    1,
    4
  ),
  top = TRUE
)

#----------------------------------------------------------
# Set ordinary header labels
#----------------------------------------------------------

ft <- set_header_labels(
  ft,
  Sex = "",
  Trait = "Trait",
  `S ± SE` = "",
  P_S = "P",
  `C ± SE` = "",
  P_C = "P",
  spacer = "",
  `β ± SE` = "",
  P_β = "P",
  `γ ± SE` = "",
  P_γ = "P"
)

#----------------------------------------------------------
# Create true subscripts in the second header row
#----------------------------------------------------------

ft <- compose(
  ft,
  i = 2,
  j = "S ± SE",
  part = "header",
  value = as_paragraph(
    as_chunk("S"),
    as_sub("i"),
    as_chunk(" ± SE")
  )
)

ft <- compose(
  ft,
  i = 2,
  j = "C ± SE",
  part = "header",
  value = as_paragraph(
    as_chunk("C"),
    as_sub("ii"),
    as_chunk(" ± SE")
  )
)

ft <- compose(
  ft,
  i = 2,
  j = "β ± SE",
  part = "header",
  value = as_paragraph(
    as_chunk("β"),
    as_sub("i"),
    as_chunk(" ± SE")
  )
)

ft <- compose(
  ft,
  i = 2,
  j = "γ ± SE",
  part = "header",
  value = as_paragraph(
    as_chunk("γ ± SE")
  )
)

# Italicize P in the column headings
for (p_col in c("P_S", "P_C", "P_β", "P_γ")) {
  
  ft <- compose(
    ft,
    i = 2,
    j = p_col,
    part = "header",
    value = as_paragraph(
      as_i("P")
    )
  )
}

#----------------------------------------------------------
# Remove default borders
#----------------------------------------------------------

ft <- border_remove(ft)

# Columns within the differential and gradient sections
diff_cols <- which(
  cols %in% c(
    "S ± SE",
    "P_S",
    "C ± SE",
    "P_C"
  )
)

grad_cols <- which(
  cols %in% c(
    "β ± SE",
    "P_β",
    "γ ± SE",
    "P_γ"
  )
)

# Explicitly exclude spacer column
diff_cols <- setdiff(diff_cols, spacer_col)
grad_cols <- setdiff(grad_cols, spacer_col)

#----------------------------------------------------------
# Header borders
#----------------------------------------------------------

# Top rule
ft <- border(
  ft,
  i = 1,
  j = seq_along(cols),
  border.top = fp_border(width = 1.5),
  part = "header"
)

# Underline spanning headings separately
ft <- border(
  ft,
  i = 1,
  j = diff_cols,
  border.bottom = fp_border(width = 1),
  part = "header"
)

ft <- border(
  ft,
  i = 1,
  j = grad_cols,
  border.bottom = fp_border(width = 1),
  part = "header"
)

# Full rule below the second header row
ft <- border(
  ft,
  i = 2,
  j = seq_along(cols),
  border.bottom = fp_border(width = 1.5),
  part = "header"
)

# Remove borders around the spacer column
ft <- border(
  ft,
  j = spacer_col,
  border.left = fp_border(width = 0),
  border.right = fp_border(width = 0),
  part = "all"
)

#----------------------------------------------------------
# Merge section-heading rows
#----------------------------------------------------------

section_rows <- which(
  table_combined$Sex %in% c(
    "(a) Females",
    "(b) Males"
  )
)

for (i in section_rows) {
  
  ft <- merge_at(
    ft,
    i = i,
    j = seq_along(cols),
    part = "body"
  )
}

#----------------------------------------------------------
# Alignment
#----------------------------------------------------------

ft <- ft %>%
  align(
    part = "header",
    align = "center"
  ) %>%
  align(
    j = c("Sex", "Trait"),
    align = "left",
    part = "all"
  ) %>%
  align(
    j = setdiff(seq_along(cols), c(1, 2, spacer_col)),
    align = "center",
    part = "all"
  )

# Left-align section headers
ft <- align(
  ft,
  i = section_rows,
  j = seq_along(cols),
  align = "left",
  part = "body"
)

# Indent trait names
trait_rows <- setdiff(
  seq_len(nrow(table_combined)),
  section_rows
)

ft <- padding(
  ft,
  i = trait_rows,
  j = "Trait",
  padding.left = 10,
  part = "body"
)

#----------------------------------------------------------
# Fonts and emphasis
#----------------------------------------------------------

ft <- ft %>%
  bold(part = "header") %>%
  bold(
    i = section_rows,
    bold = TRUE,
    part = "body"
  ) %>%
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  font(
    fontname = "Times New Roman",
    part = "all"
  )

#----------------------------------------------------------
# Bold statistically significant results
#----------------------------------------------------------

alpha <- 0.05

sig_S <- which(
  !is.na(table_combined$P_S_num) &
    table_combined$P_S_num < alpha
)

sig_C <- which(
  !is.na(table_combined$P_C_num) &
    table_combined$P_C_num < alpha
)

sig_beta <- which(
  !is.na(table_combined$P_β_num) &
    table_combined$P_β_num < alpha
)

sig_gamma <- which(
  !is.na(table_combined$P_γ_num) &
    table_combined$P_γ_num < alpha
)

if (length(sig_S) > 0) {
  
  ft <- bold(
    ft,
    i = sig_S,
    j = c("S ± SE", "P_S"),
    part = "body"
  )
}

if (length(sig_C) > 0) {
  
  ft <- bold(
    ft,
    i = sig_C,
    j = c("C ± SE", "P_C"),
    part = "body"
  )
}

if (length(sig_beta) > 0) {
  
  ft <- bold(
    ft,
    i = sig_beta,
    j = c("β ± SE", "P_β"),
    part = "body"
  )
}

if (length(sig_gamma) > 0) {
  
  ft <- bold(
    ft,
    i = sig_gamma,
    j = c("γ ± SE", "P_γ"),
    part = "body"
  )
}

#----------------------------------------------------------
# Final sizing and borders
#----------------------------------------------------------

ft <- autofit(ft)

# Narrow spacer between differential and gradient sections
ft <- width(
  ft,
  j = spacer_col,
  width = 0.05
)

# Bottom rule
ft <- border(
  ft,
  i = nrow(table_combined),
  j = seq_along(cols),
  border.bottom = fp_border(width = 1.5),
  part = "body"
)

# Reduce cell padding
ft <- padding(
  ft,
  padding.top = 2,
  padding.bottom = 2,
  part = "all"
)

# Display the table
ft

# =========================
# female fecundity
# =========================

females_clean <- females |>
  filter(!is.na(eggs), total_body > 0)

# Fit negative binomial GLM to model egg load vs body length
female.fecundity <- glm.nb(eggs ~ log(total_body), data=females_clean)
summary(female.fecundity) # strong fecundity selection

newdat <- data.frame(
  total_body = seq(min(females_clean$total_body),
                 max(females_clean$total_body),
                 length.out = 100)
)

pred <- predict(female.fecundity, newdata = newdat,
                type = "link", se.fit = TRUE)

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
#       labs(x = "Body length (mm, log scale)",
#                y = "Egg load (number of eggs)") +
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
# ggsave("figure_2.png", width = 8, height = 8, dpi = 600)


# =========================
# assortative mating
# =========================

# Restructure data to compare male and female traits in pairs
pairs <- morph_data %>%
  filter(!is.na(pair), pair != "NA") %>%
  pivot_wider(
    id_cols = pair,
    names_from = sex,
    values_from = c(total_body, total_leg, mated)
  ) %>%
  drop_na(total_body_f, total_body_m)

stopifnot(nrow(pairs) == 25)

# Test for assortative mating using standardized linear regression
assort <- lm(scale(total_body_f) ~ scale(total_body_m), data = pairs)
summary(assort)

# The standardized regression slope is equal to the correlation
obs_slope <- cor(pairs$total_body_m, pairs$total_body_f)

# Two-sided randomization test using slopes on the same standardized scale
set.seed(2029)

null_slopes <- replicate(
  10000,
  cor(
    pairs$total_body_m,
    sample(pairs$total_body_f, replace = FALSE)
  )
)

p_val <- (
  sum(abs(null_slopes) >= abs(obs_slope)) + 1
) / (
  length(null_slopes) + 1
)

obs_slope
p_val

# ============================================================
# Sex-specific comparison of multivariate selection gradients
# ============================================================

# Required packages
library(dplyr)
library(tibble)


# ------------------------------------------------------------
# 1. Check and prepare the data
# ------------------------------------------------------------

# Check that sex contains only the expected values
sex_values <- unique(tolower(trimws(morph_data$sex)))

if (!all(sex_values %in% c("f", "m"))) {
  stop(
    "Unexpected values found in sex: ",
    paste(setdiff(sex_values, c("f", "m")), collapse = ", ")
  )
}


selection_data <- morph_data %>%
  transmute(
    id = as.character(id),
    pair = pair,
    
    sex = case_when(
      tolower(trimws(sex)) == "f" ~ "Female",
      tolower(trimws(sex)) == "m" ~ "Male",
      TRUE ~ NA_character_
    ),
    
    mated = as.integer(mated),
    
    total_body = as.numeric(total_body),
    total_leg  = as.numeric(total_leg)
  ) %>%
  
  # Retain complete cases for variables used in the analysis
  filter(
    !is.na(id),
    !is.na(sex),
    !is.na(mated),
    !is.na(total_body),
    !is.na(total_leg)
  ) %>%
  
  mutate(
    sex = factor(
      sex,
      levels = c("Female", "Male")
    )
  )


# Confirm that mating success is coded 0/1
if (!all(selection_data$mated %in% c(0L, 1L))) {
  stop("The variable 'mated' must contain only 0 and 1.")
}


# ------------------------------------------------------------
# 2. Create clusters for paired observations
# ------------------------------------------------------------

# A paired male and female have the same pair number.
# Unpaired individuals are each assigned their own cluster.
selection_data <- selection_data %>%
  mutate(
    cluster_id = if_else(
      !is.na(pair),
      paste0("pair_", pair),
      paste0("unpaired_", id)
    )
  )


# ------------------------------------------------------------
# 3. Calculate sex-specific means and standard deviations
# ------------------------------------------------------------

scaling_values <- selection_data %>%
  group_by(sex) %>%
  summarise(
    body_mean = mean(total_body),
    body_sd   = sd(total_body),
    
    leg_mean = mean(total_leg),
    leg_sd   = sd(total_leg),
    
    .groups = "drop"
  )

scaling_values


# Check that there is variation within each sex
if (
  any(scaling_values$body_sd <= 0) ||
  any(scaling_values$leg_sd <= 0)
) {
  stop("Body or leg length has no within-sex variation.")
}


# ------------------------------------------------------------
# 4. Standardize morphology within sex
# ------------------------------------------------------------

selection_data <- selection_data %>%
  left_join(
    scaling_values,
    by = "sex"
  ) %>%
  mutate(
    # Within-sex z scores
    body_z = (total_body - body_mean) / body_sd,
    leg_z  = (total_leg - leg_mean) / leg_sd,
    
    # Second-order terms
    body_z2 = body_z^2,
    leg_z2  = leg_z^2,
    
    # Correlational term
    body_leg = body_z * leg_z
  )


# Check the transformation
selection_data %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    body_z_mean = mean(body_z),
    body_z_sd = sd(body_z),
    leg_z_mean = mean(leg_z),
    leg_z_sd = sd(leg_z),
    .groups = "drop"
  )


# Sample sizes and pairing rates
sample_summary <- selection_data %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    n_mated = sum(mated),
    n_unmated = sum(mated == 0),
    mating_rate = mean(mated),
    .groups = "drop"
  )

sample_summary


# ------------------------------------------------------------
# 5. Fit the pooled selection-surface models
# ------------------------------------------------------------

# Common surface:
# morphology has the same relationship with pairing success
# in males and females.
common_surface_model <- glm(
  mated ~
    sex +
    body_z +
    leg_z +
    body_z2 +
    leg_z2 +
    body_leg,
  family = binomial(link = "logit"),
  data = selection_data,
  control = glm.control(maxit = 100)
)


# Full sex-specific surface:
# sex interacts with every linear, quadratic, and
# correlational term.
sex_specific_model <- glm(
  mated ~
    sex * (
      body_z +
        leg_z +
        body_z2 +
        leg_z2 +
        body_leg
    ),
  family = binomial(link = "logit"),
  data = selection_data,
  control = glm.control(maxit = 100)
)


# Reduced model in which the correlational gradient is
# constrained to be the same in males and females.
#
# All other sex-specific linear and quadratic effects remain.
no_sex_correlational_model <- glm(
  mated ~
    sex * (
      body_z +
        leg_z +
        body_z2 +
        leg_z2
    ) +
    body_leg,
  family = binomial(link = "logit"),
  data = selection_data,
  control = glm.control(maxit = 100)
)


summary(sex_specific_model)


# ------------------------------------------------------------
# 6. Basic model diagnostics
# ------------------------------------------------------------

model_diagnostics <- tibble(
  converged = sex_specific_model$converged,
  iterations = sex_specific_model$iter,
  maximum_absolute_coefficient =
    max(abs(coef(sex_specific_model))),
  minimum_fitted_probability =
    min(fitted(sex_specific_model)),
  maximum_fitted_probability =
    max(fitted(sex_specific_model))
)

model_diagnostics


if (!sex_specific_model$converged) {
  warning("The full sex-specific model did not converge.")
}

if (max(abs(coef(sex_specific_model))) > 10) {
  warning(
    "One or more coefficients are very large. ",
    "Inspect the model for sparse data or separation."
  )
}


# ------------------------------------------------------------
# 7. Model-based likelihood-ratio tests
# ------------------------------------------------------------

# Omnibus test:
# Does any component of the selection surface differ by sex?
overall_sex_LRT <- anova(
  common_surface_model,
  sex_specific_model,
  test = "Chisq"
)

overall_sex_LRT


# Specific test:
# Does the body-by-leg correlational term differ by sex?
correlational_sex_LRT <- anova(
  no_sex_correlational_model,
  sex_specific_model,
  test = "Chisq"
)

correlational_sex_LRT


# ------------------------------------------------------------
# 8. Pair-clustered covariance matrix
# ------------------------------------------------------------

# Paired males and females share a mating event and therefore
# should not necessarily be treated as fully independent.
#
# Unpaired individuals are singleton clusters.
clustered_glm_vcov <- function(model, cluster) {
  design_matrix <- model.matrix(model)
  response_residual <- residuals(model, type = "response")
  observation_scores <- design_matrix * response_residual
  cluster_scores <- rowsum(
    observation_scores,
    group = cluster,
    reorder = FALSE
  )

  information_matrix <- crossprod(
    design_matrix,
    design_matrix * model$weights
  )
  bread <- solve(information_matrix)
  meat <- crossprod(cluster_scores)

  n_observations <- nrow(design_matrix)
  n_coefficients <- ncol(design_matrix)
  n_clusters <- nrow(cluster_scores)
  finite_sample_correction <-
    (n_clusters / (n_clusters - 1)) *
    ((n_observations - 1) / (n_observations - n_coefficients))

  finite_sample_correction * bread %*% meat %*% bread
}

clustered_vcov <- clustered_glm_vcov(
  sex_specific_model,
  selection_data$cluster_id
)

clustered_standard_errors <- sqrt(diag(clustered_vcov))
clustered_statistics <- coef(sex_specific_model) / clustered_standard_errors

clustered_coefficients <- tibble(
  term = names(coef(sex_specific_model)),
  estimate = unname(coef(sex_specific_model)),
  std.error = unname(clustered_standard_errors),
  statistic = unname(clustered_statistics),
  p.value = 2 * pnorm(abs(clustered_statistics), lower.tail = FALSE)
)

clustered_coefficients


# ------------------------------------------------------------
# 9. Omnibus cluster-robust test of all sex interactions
# ------------------------------------------------------------

coefficient_names <- names(coef(sex_specific_model))

sex_interaction_terms <- coefficient_names[
  grepl("sexMale", coefficient_names, fixed = TRUE) &
    grepl(":", coefficient_names, fixed = TRUE)
]

sex_interaction_terms


# Construct a hypothesis matrix testing all sex interactions
# simultaneously.
K_all_sex_interactions <- matrix(
  0,
  nrow = length(sex_interaction_terms),
  ncol = length(coefficient_names),
  dimnames = list(
    sex_interaction_terms,
    coefficient_names
  )
)

for (i in seq_along(sex_interaction_terms)) {
  K_all_sex_interactions[
    i,
    sex_interaction_terms[i]
  ] <- 1
}


interaction_estimates <- drop(
  K_all_sex_interactions %*% coef(sex_specific_model)
)
interaction_covariance <-
  K_all_sex_interactions %*%
  clustered_vcov %*%
  t(K_all_sex_interactions)
overall_clustered_chisq <- drop(
  t(interaction_estimates) %*%
  solve(interaction_covariance, interaction_estimates)
)
overall_clustered_df <- qr(interaction_covariance)$rank
overall_clustered_test <- tibble(
  statistic = overall_clustered_chisq,
  df = overall_clustered_df,
  p.value = pchisq(
    overall_clustered_chisq,
    df = overall_clustered_df,
    lower.tail = FALSE
  )
)

overall_clustered_test


# ------------------------------------------------------------
# 10. Locate interaction coefficient names safely
# ------------------------------------------------------------

# Depending on how R constructs the model matrix, an interaction
# can be named sexMale:body_leg or body_leg:sexMale.
interaction_name <- function(variable) {
  
  possible_names <- c(
    paste0("sexMale:", variable),
    paste0(variable, ":sexMale")
  )
  
  found_name <- possible_names[
    possible_names %in% names(coef(sex_specific_model))
  ]
  
  if (length(found_name) != 1) {
    stop(
      "Could not uniquely identify the sex interaction for ",
      variable
    )
  }
  
  found_name
}


interaction_name("body_z")
interaction_name("leg_z")
interaction_name("body_z2")
interaction_name("leg_z2")
interaction_name("body_leg")


# ------------------------------------------------------------
# 11. Direct robust test of the correlational-gradient difference
# ------------------------------------------------------------

correlational_interaction_term <- interaction_name("body_leg")

correlational_interaction_test <- clustered_coefficients %>%
  filter(term == correlational_interaction_term)

correlational_interaction_test

# ------------------------------------------------------------
# 12. Function for linear combinations of coefficients
# ------------------------------------------------------------

estimate_linear_combination <- function(
    weights,
    multiplier = 1,
    model = sex_specific_model,
    vcov_matrix = clustered_vcov,
    confidence_level = 0.95
) {
  
  model_coefficients <- coef(model)
  coefficient_names <- names(model_coefficients)
  
  missing_terms <- setdiff(
    names(weights),
    coefficient_names
  )
  
  if (length(missing_terms) > 0) {
    stop(
      "Terms not found in the model: ",
      paste(missing_terms, collapse = ", ")
    )
  }
  
  contrast_vector <- setNames(
    rep(0, length(model_coefficients)),
    coefficient_names
  )
  
  contrast_vector[names(weights)] <- weights
  
  raw_estimate <- sum(
    contrast_vector * model_coefficients
  )
  
  raw_variance <- drop(
    t(contrast_vector) %*%
      vcov_matrix %*%
      contrast_vector
  )
  
  # Protect against a tiny negative value caused by
  # floating-point precision.
  raw_variance <- max(raw_variance, 0)
  
  raw_se <- sqrt(raw_variance)
  
  estimate <- multiplier * raw_estimate
  std_error <- abs(multiplier) * raw_se
  
  z_value <- estimate / std_error
  
  p_value <- 2 * pnorm(
    abs(z_value),
    lower.tail = FALSE
  )
  
  alpha <- 1 - confidence_level
  critical_value <- qnorm(1 - alpha / 2)
  
  tibble(
    estimate = estimate,
    std.error = std_error,
    conf.low = estimate - critical_value * std_error,
    conf.high = estimate + critical_value * std_error,
    z = z_value,
    p.value = p_value
  )
}


# ------------------------------------------------------------
# 13. Function for a sex-specific gradient
# ------------------------------------------------------------

estimate_sex_gradient <- function(
    sex_name,
    model_term,
    gradient_name,
    multiplier = 1
) {
  
  if (sex_name == "Female") {
    
    weights <- setNames(
      1,
      model_term
    )
    
  } else if (sex_name == "Male") {
    
    weights <- setNames(
      c(1, 1),
      c(
        model_term,
        interaction_name(model_term)
      )
    )
    
  } else {
    
    stop("sex_name must be 'Female' or 'Male'.")
  }
  
  estimate_linear_combination(
    weights = weights,
    multiplier = multiplier
  ) %>%
    mutate(
      sex = sex_name,
      gradient = gradient_name,
      .before = 1
    )
}


# ------------------------------------------------------------
# 14. Sex-specific gradient table
# ------------------------------------------------------------

sex_specific_gradients <- bind_rows(
  
  # Female gradients
  estimate_sex_gradient(
    "Female",
    "body_z",
    "beta_body"
  ),
  
  estimate_sex_gradient(
    "Female",
    "leg_z",
    "beta_leg"
  ),
  
  estimate_sex_gradient(
    "Female",
    "body_z2",
    "gamma_body",
    multiplier = 2
  ),
  
  estimate_sex_gradient(
    "Female",
    "leg_z2",
    "gamma_leg",
    multiplier = 2
  ),
  
  estimate_sex_gradient(
    "Female",
    "body_leg",
    "gamma_body_leg"
  ),
  
  # Male gradients
  estimate_sex_gradient(
    "Male",
    "body_z",
    "beta_body"
  ),
  
  estimate_sex_gradient(
    "Male",
    "leg_z",
    "beta_leg"
  ),
  
  estimate_sex_gradient(
    "Male",
    "body_z2",
    "gamma_body",
    multiplier = 2
  ),
  
  estimate_sex_gradient(
    "Male",
    "leg_z2",
    "gamma_leg",
    multiplier = 2
  ),
  
  estimate_sex_gradient(
    "Male",
    "body_leg",
    "gamma_body_leg"
  )
) %>%
  dplyr::select(
    sex,
    gradient,
    estimate,
    std.error,
    conf.low,
    conf.high,
    z,
    p.value
  )


sex_specific_gradients

# ------------------------------------------------------------
# 15. Function for male-minus-female differences
# ------------------------------------------------------------

estimate_sex_difference <- function(
    model_term,
    gradient_name,
    multiplier = 1
) {
  
  weights <- setNames(
    1,
    interaction_name(model_term)
  )
  
  estimate_linear_combination(
    weights = weights,
    multiplier = multiplier
  ) %>%
    mutate(
      contrast = "Male - Female",
      gradient = gradient_name,
      .before = 1
    )
}


# ------------------------------------------------------------
# 16. Male-minus-female gradient differences
# ------------------------------------------------------------

sex_gradient_differences <- bind_rows(
  
  estimate_sex_difference(
    "body_z",
    "beta_body"
  ),
  
  estimate_sex_difference(
    "leg_z",
    "beta_leg"
  ),
  
  estimate_sex_difference(
    "body_z2",
    "gamma_body",
    multiplier = 2
  ),
  
  estimate_sex_difference(
    "leg_z2",
    "gamma_leg",
    multiplier = 2
  ),
  
  estimate_sex_difference(
    "body_leg",
    "gamma_body_leg"
  )
) %>%
  dplyr::select(
    contrast,
    gradient,
    estimate,
    std.error,
    conf.low,
    conf.high,
    z,
    p.value
  )


sex_gradient_differences

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

# Load male-male combat data
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

# Calculate size differences between winners and losers
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

# Model the effect of size differences on winning probability
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



