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
library(brms)
library(broom.mixed)
library(parameters)
library(posterior)
library(gsg)
library(mgcv)
library(MASS)

#install.packages("~/Documents/Action/gsg_2.0.tar", repos = NULL, type = "source")

#### read in data file ####
morph_data <- read.csv(file="data/raw/weevils.csv", header=TRUE, sep=",", dec=".") %>%
  as_tibble() %>%
  print(n=200)

morph_data <- morph_data %>%
  mutate(fem=rowMeans(dplyr::select(., l_fem, r_fem), na.rm = TRUE)) %>%
  mutate(tib=rowMeans(dplyr::select(., l_tib, r_tib), na.rm = TRUE)) %>%
  mutate(total_body=tot_abdo+thorax) %>%
  mutate(total_leg=fem+tib)

# =============================
# sexual size dimorphism
# =============================
# make sure sex is a factor
morph_data$sex <- factor(morph_data$sex)

# total body length
model_body <- lm(total_body ~ sex, data = morph_data)
summary(model_body)

# total leg length
model_leg <- lm(total_leg ~ sex, data = morph_data)
summary(model_leg)

manova_model <- manova(cbind(total_body, total_leg) ~ sex, data = morph_data)
summary(manova_model)

ggplot(morph_data, aes(x = sex, y = total_body)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)

ggplot(morph_data, aes(x = sex, y = total_leg)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5)

# =============================
# selection analysis
# =============================

# male selection analysis
males <- morph_data %>%
  dplyr::filter(sex=="m")

male.gam<-gam(mated~s(total_body)+s(total_leg),family=binomial(link="logit"),method="GCV.Cp", data=males)
male.grad<-gam.gradients(mod = male.gam,phenotype=c("total_body","total_leg"),se.method="boot.para",n.boot=1000,standardize=TRUE)

male.estimates<-data.frame(male.grad$ests)

## need to double the (quadratic) gamma gradient and SEs
male.estimates[3,1]<-male.estimates[3,1]*2
male.estimates[4,1]<-male.estimates[4,1]*2
male.estimates[3,2]<-male.estimates[3,2]*2
male.estimates[4,2]<-male.estimates[4,2]*2
save(male.estimates, file = "male.estimates.RData")

#male selection differentials
m.body <- moments.differentials(z=males$total_body,W=males$mated, n.boot=10000, standardized=TRUE)
m.fem <- moments.differentials(z=males$fem,W=males$mated, n.boot=10000, standardized=TRUE)
m.tib <- moments.differentials(z=males$tib,W=males$mated, n.boot=10000, standardized=TRUE)

# plot body length
m.body <- fitness.landscape(mod=male.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(m.body$points[,1],m.body$Wbar,type="l", ylim=c(0,1),xlab="Body length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(m.body$points[,1],m.body$WbarPI[1,],lty=2)
lines(m.body$points[,1],m.body$WbarPI[2,],lty=2)

# plot leg length
m.leg <- fitness.landscape(mod=male.gam, phenotype="total_leg", covariates=c("total_body"), PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(m.leg$points[,1],m.leg$Wbar,type="l", ylim=c(0,1),xlab="Hind leg length (mm)", ylab="Mating success", cex.lab=2, cex.axis=1.5)
lines(m.leg$points[,1],m.leg$WbarPI[1,],lty=2)
lines(m.leg$points[,1],m.leg$WbarPI[2,],lty=2)

dev.off()
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0), # move plot to the right and up
    mgp = c(2, 1, 0) # move axis labels closer to axis
) 
vis.gam(male.gam,view=c("total_body","total_leg"),color="heat",plot.type = "persp",  type="response",theta=155,
        xlab="\nBody length (mm)", ylab="\nLeg length. (mm)",zlab="\nFitness")
male.proj<-recordPlot()

# female selection analysis
females <- morph_data %>%
  dplyr::filter(sex=="f")

female.gam<-gam(mated~s(total_body)+s(total_leg),family=binomial(link="logit"),method="GCV.Cp", data=females)
female.grad<-gam.gradients(mod = female.gam, phenotype=c("total_body","total_leg"),se.method="boot.para",n.boot=1000,standardize=TRUE)

female.estimates<-data.frame(female.grad$ests)

## need to double the (quadratic) gamma gradient and SEs
female.estimates[3,1]<-female.estimates[3,1]*2
female.estimates[4,1]<-female.estimates[4,1]*2
female.estimates[3,2]<-female.estimates[3,2]*2
female.estimates[4,2]<-female.estimates[4,2]*2
save(female.estimates, file = "female.estimates.RData")

# female selection differentials
f.body <- moments.differentials(z=females$total_body,W=females$mated, n.boot=10000, standardized=TRUE)
f.leg <- moments.differentials(z=females$total_leg,W=females$mated, n.boot=10000, standardized=TRUE)
-0.57^2
# plot abdoment length
f.body <- fitness.landscape(mod=female.gam,phenotype="total_body",covariates=c("total_leg"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.body$points[,1],f.body$Wbar,type="l", ylim=c(0,1),xlab="Body length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.body$points[,1],f.body$WbarPI[1,],lty=2)
lines(f.body$points[,1],f.body$WbarPI[2,],lty=2)

# plot femur length
f.leg <- fitness.landscape(mod=female.gam,phenotype="total_leg",covariates=c("total_body"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.leg$points[,1],f.leg$Wbar,type="l", ylim=c(0,1),xlab="Hind leg length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.leg$points[,1],f.leg$WbarPI[1,],lty=2)
lines(f.leg$points[,1],f.leg$WbarPI[2,],lty=2)

dev.off()

# =========================
# female fecundity
# =========================

female.fecundity <- glm.nb(eggs ~ log(total_body), data=females)
summary(female.fecundity) # strong fecundity selection

newdat <- data.frame(
  total_body = seq(min(females$total_body),
                 max(females$total_body),
                 length.out = 100)
)

pred <- predict(female.fecundity, newdata = newdat,
                type = "link", se.fit = TRUE)

newdat$fit <- exp(pred$fit)
newdat$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
newdat$upr <- exp(pred$fit + 1.96 * pred$se.fit)

female.fecundity.plot <- ggplot(females, aes(x = total_body, y = eggs)) +
      geom_point(alpha = 0.5) +
      geom_line(data = newdat, aes(x = total_body, y = fit),
            color = "blue", linewidth = 1) +
      geom_ribbon(data = newdat,
              aes(x = total_body, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              alpha = 0.2) +
      scale_x_log10() +
      labs(x = "Body size (mm, log scale)",
               y = "Fecundity (number of eggs)") +
      theme_classic()

# =========================
# assortative mating
# =========================

pairs <- morph_data %>%
  filter(pair!="NA")%>%
  pivot_wider(
    id_cols = pair,
    names_from = sex,
    values_from = c(fem, tib, tot_abdo, thorax, mated)
  ) %>%
  rename(
    female_fem = fem_f,
    male_fem   = fem_m,
    female_tib = tib_f,
    male_tib   = tib_m
  )

assort <- lm(scale(tot_abdo_f) ~ scale(tot_abdo_m), data = pairs)
summary(assort)

obs_slope <- coef(assort)[2]

null_slopes <- replicate(10000, {
  perm_females <- sample(pairs$tot_abdo_f)
  coef(lm(perm_females ~ pairs$tot_abdo_m))[2]
})

p_val <- mean(abs(null_slopes) >= abs(obs_slope)) # no assortative mating

# =========================
# male mate choice
# =========================

male_mate_choice <- read.csv(file="data/raw/male_mate_choice.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

choice <- male_mate_choice %>%
  left_join(
    females %>% dplyr::select(id, total_body, total_leg),
    by = c("chosen_female_id" = "id")
  ) %>%
  rename(
    total_body_chosen = total_body,
    leg_chosen        = total_leg
  ) %>%
  left_join(
    females %>% dplyr::select(id, total_body, total_leg),
    by = c("unchosen_fem_id" = "id")
  ) %>%
  rename(
    total_body_unchosen = total_body,
    leg_unchosen        = total_leg
  )

choice$chose_larger <- ifelse(choice$total_body_chosen > choice$total_body_unchosen, 1, 0)
choice$size_diff <- choice$total_body_chosen - choice$total_body_unchosen

table(choice$chose_larger)

binom.test(sum(choice$chose_larger), nrow(choice), p = 0.5)

# =========================
# male-male competition
# =========================

male_male <- read.csv(file="data/raw/male_male_combat.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

# dyadic
combat2 <- male_male %>%
  left_join(
    males %>% dplyr::select(id, total_body, fem),
    by = c("winner_id" = "id")
  ) %>%
  rename(
    total_body_win = total_body,
    fem_win        = fem
  ) %>%
  left_join(
    males %>% dplyr::select(id, total_body, fem),
    by = c("loser_id" = "id")
  ) %>%
  rename(
    total_body_lose = total_body,
    fem_lose        = fem
  )

combat2$won_larger <- ifelse(combat2$total_body_win > combat2$total_body_lose, 1, 0)
combat2$size_diff <- combat2$total_body_win - combat2$total_body_lose

binom.test(sum(combat2$won_larger), nrow(combat2), p = 0.5)

# combat2 <- combat2 %>%
#   mutate(
#     body_diff = total_body_win - total_body_lose,
#     fem_diff  = fem_win - fem_lose
#   ) # positive = winner is bigger
# 
# df <- combat2
# 
# set.seed(123)
# 
# flip <- rbinom(nrow(df), 1, 0.5)
# 
# df$id_i <- ifelse(flip == 1, df$winner_id, df$loser_id)
# df$id_j <- ifelse(flip == 1, df$loser_id, df$winner_id)
# 
# df$fem_i <- ifelse(flip == 1, df$fem_win, df$fem_lose)
# df$fem_j <- ifelse(flip == 1, df$fem_lose, df$fem_win)
# 
# df$body_i <- ifelse(flip == 1, df$total_body_win, df$total_body_lose)
# df$body_j <- ifelse(flip == 1, df$total_body_lose, df$total_body_win)
# 
# df$win <- ifelse(df$id_i == df$winner_id, 1, 0)
# 
# df$diff_fem  <- df$fem_i - df$fem_j
# df$diff_body <- df$body_i - df$body_j
# 
# model.dyad <- glm(win ~ scale(diff_fem) + scale(diff_body),
#              family = binomial, data = df)
# 
# summary(model.dyad)
# 
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
# winners <- data.frame(
#   id   = combat2$winner_id,
#   fem  = combat2$fem_win,
#   body = combat2$total_body_win,
#   win  = 1
# )
# 
# losers <- data.frame(
#   id   = combat2$loser_id,
#   fem  = combat2$fem_lose,
#   body = combat2$total_body_lose,
#   win  = 0
# )
# 
# indiv <- rbind(winners, losers)
# 
# indiv$fem_z  <- scale(indiv$fem)
# indiv$body_z <- scale(indiv$body)
# 
# model_sel <- glm(win ~ fem_z + body_z, family = binomial, data=indiv)
# summary(model_sel)







## ---- end