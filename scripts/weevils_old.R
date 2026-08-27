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
library(broman)
library(scales)
# install.packages("~/Documents/Action/gsg_2.0.tar", repos = NULL, type = "source")

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

morph_data %>%
  filter(sex=="m") %>%
  filter(total_leg<6)

# Calculate summary statistics by sex
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
# Fit linear models to test for sex differences in body and leg length
model_body <- lm(total_body ~ sex, data = morph_data)
summary(model_body)

# total leg length
model_leg <- lm(total_leg ~ sex, data = morph_data)
summary(model_leg)

# Test for multivariate sexual dimorphism using MANOVA
manova_model <- manova(cbind(total_body, total_leg) ~ sex, data = morph_data)
summary(manova_model)

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

male.grad <- gam.gradients(
  mod = male.gam,
  phenotype = c("total_body", "total_leg"),
  se.method = "boot.para",
  n.boot = 10000,
  standardized = TRUE,
  refit.smooth = FALSE
)

female.grad <- gam.gradients(
  mod = female.gam,
  phenotype = c("total_body", "total_leg"),
  se.method = "boot.para",
  n.boot = 10000,
  standardized = TRUE,
  refit.smooth = FALSE
)

male.estimates <- as.data.frame(male.grad$ests)
female.estimates <- as.data.frame(female.grad$ests)

saveRDS(male.estimates, "data/processed/male.estimates.rds")
saveRDS(female.estimates, "data/processed/female.estimates.rds")

# Figure 3: explicit 95% parametric-bootstrap intervals
set.seed(2026)

m.body.fitness <- fitness.landscape(
  mod = male.gam,
  phenotype = "total_body",
  covariates = "total_leg",
  PI.method = "boot.para",
  PI.interval = c(0.025, 0.975),
  n.boot = 1000,
  refit.smooth = FALSE
)

m.leg.fitness <- fitness.landscape(
  mod = male.gam,
  phenotype = "total_leg",
  covariates = "total_body",
  PI.method = "boot.para",
  PI.interval = c(0.025, 0.975),
  n.boot = 1000,
  refit.smooth = FALSE
)

f.body.fitness <- fitness.landscape(
  mod = female.gam,
  phenotype = "total_body",
  covariates = "total_leg",
  PI.method = "boot.para",
  PI.interval = c(0.025, 0.975),
  n.boot = 1000,
  refit.smooth = FALSE
)

f.leg.fitness <- fitness.landscape(
  mod = female.gam,
  phenotype = "total_leg",
  covariates = "total_body",
  PI.method = "boot.para",
  PI.interval = c(0.025, 0.975),
  n.boot = 1000,
  refit.smooth = FALSE
)

saveRDS(m.body.fitness, "data/processed/m.body.fitness.rds")
saveRDS(m.leg.fitness,  "data/processed/m.leg.fitness.rds")
saveRDS(f.body.fitness, "data/processed/f.body.fitness.rds")
saveRDS(f.leg.fitness,  "data/processed/f.leg.fitness.rds")

#----------------------------------------------------------
# 1. Extract the observations used in the model
#----------------------------------------------------------

male_data <- model.frame(male.gam)

male_data <- male_data[
  complete.cases(
    male_data$total_body,
    male_data$total_leg
  ),
]

#----------------------------------------------------------
# 2. Construct a prediction grid
#----------------------------------------------------------

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
  total_leg  = y_seq
)

pred_grid$probability <- predict(
  male.gam,
  newdata = pred_grid,
  type = "response"
)

#----------------------------------------------------------
# 3. Identify the convex hull of the observations
#----------------------------------------------------------

hull_id <- chull(
  male_data$total_body,
  male_data$total_leg
)

hull_x <- male_data$total_body[hull_id]
hull_y <- male_data$total_leg[hull_id]

# Close the polygon
hull_x <- c(hull_x, hull_x[1])
hull_y <- c(hull_y, hull_y[1])

# Determine which prediction cells fall within the hull
inside_hull <- sp::point.in.polygon(
  point.x = pred_grid$total_body,
  point.y = pred_grid$total_leg,
  pol.x   = hull_x,
  pol.y   = hull_y
) > 0

# Mask unsupported trait combinations
pred_grid$probability[!inside_hull] <- NA_real_

# Convert predictions to the matrix required by image()
z_matrix <- matrix(
  pred_grid$probability,
  nrow = length(x_seq),
  ncol = length(y_seq)
)

#----------------------------------------------------------
# 4. Draw the figure
#----------------------------------------------------------

png(
  "figure_4.png",
  width = 2200,
  height = 1900,
  res = 300,
  bg = "white"
)

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
  ylab = "Leg length (mm)",
  main = ""
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

points(
  male_data$total_body,
  male_data$total_leg,
  pch = 16,
  cex = 0.4,
  col = adjustcolor("black", alpha.f = 0.4)
)

box()

dev.off()

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
       ylim=c(0,0.8),
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
        trait == "total_leg"  ~ "leg length",
        TRUE                  ~ trait
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

#----------------------------------------------------------
# Format estimates and standard errors
#----------------------------------------------------------

final_table <- bind_rows(
  female_table,
  male_table
) %>%
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
    as_chunk("γ"),
    as_sub("ii"),
    as_chunk(" ± SE")
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

str(female.fecundity)

newdat <- data.frame(
  total_body = seq(min(females$total_body),
                 max(females$total_body),
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
  filter(pair!="NA")%>%
  pivot_wider(
    id_cols = pair,
    names_from = sex,
    values_from = c(total_body, total_leg, mated)
  )

# Test for assortative mating using linear regression
assort <- lm(scale(total_body_f) ~ scale(total_body_m), data = pairs)
summary(assort)

obs_slope <- coef(assort)[2]

# Permutation test to assess significance of assortment
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







