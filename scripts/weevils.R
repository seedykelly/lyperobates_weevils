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
  as.data.frame()

male_male <- read.csv(file="data/raw/male_male_combat.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

morph_data <- morph_data %>%
  mutate(fem=rowMeans(dplyr::select(., l_fem, r_fem), na.rm = TRUE)) %>%
  mutate(tib=rowMeans(dplyr::select(., l_tib, r_tib), na.rm = TRUE)) %>%
  mutate(total_body=tot_abdo+thorax)

# =============================
# selection analysis
# =============================

# male selection analysis
males <- morph_data %>%
  dplyr::filter(sex=="m")

male.gam<-gam(mated~s(total_body)+s(fem)+s(tib),family=binomial(link="logit"),method="GCV.Cp", data=males)
male.grad<-gam.gradients(mod = male.gam,phenotype=c("total_body","fem","tib"),se.method="boot.para",n.boot=1000,standardize=TRUE)

male.estimates<-data.frame(male.grad$ests)

## need to double the (quadratic) gamma gradient
male.estimates[6,1]<-male.estimates[6,1]*2
male.estimates[7,1]<-male.estimates[7,1]*2
male.estimates[8,1]<-male.estimates[8,1]*2
male.estimates[9,1]<-male.estimates[9,1]*2
male.estimates[10,1]<-male.estimates[10,1]*2
save(male.estimates, file = "male.estimates.RData")

#male selection differentials
m.body <- moments.differentials(z=males$total_body,W=males$mated, n.boot=10000, standardized=TRUE)
m.fem <- moments.differentials(z=males$fem,W=males$mated, n.boot=10000, standardized=TRUE)
m.tib <- moments.differentials(z=males$tib,W=males$mated, n.boot=10000, standardized=TRUE)

# plot abdoment length
m.body <- fitness.landscape(mod=male.gam,phenotype="total_body",covariates=c("tib","fem"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(m.body$points[,1],m.body$Wbar,type="l", ylim=c(0,1),xlab="Body length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(m.body$points[,1],m.body$WbarPI[1,],lty=2)
lines(m.body$points[,1],m.body$WbarPI[2,],lty=2)

# plot femur length
m.femur <- fitness.landscape(mod=male.gam,phenotype="fem",covariates=c("tib","total_body"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(m.femur$points[,1],m.femur$Wbar,type="l", ylim=c(0,1),xlab="Femur length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(m.femur$points[,1],m.femur$WbarPI[1,],lty=2)
lines(m.femur$points[,1],m.femur$WbarPI[2,],lty=2)

dev.off()
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0), # move plot to the right and up
    mgp = c(2, 1, 0) # move axis labels closer to axis
) 
vis.gam(male.gam,view=c("total_body","fem"),color="heat",plot.type = "persp",  type="response",theta=155,
        xlab="\nBody length", ylab="\nFemur length",zlab="\nFitness")
male.proj<-recordPlot()

# female selection analysis
females <- morph_data %>%
  dplyr::filter(sex=="f")

female.gam<-gam(mated~s(total_body)+s(fem)+s(tib),family=binomial(link="logit"),method="GCV.Cp", data=females)
female.grad<-gam.gradients(mod = female.gam, phenotype=c("total_body","fem","tib"),se.method="boot.para",n.boot=1000,standardize=TRUE)

female.estimates<-data.frame(female.grad$ests)

## need to double the (quadratic) gamma gradient
female.estimates[6,1]<-female.estimates[6,1]*2
female.estimates[7,1]<-female.estimates[7,1]*2
female.estimates[8,1]<-female.estimates[8,1]*2
female.estimates[9,1]<-female.estimates[9,1]*2
female.estimates[10,1]<-female.estimates[10,1]*2
save(female.estimates, file = "female.estimates.RData")

# female selection differentials
f.body <- moments.differentials(z=females$tot_abdo,W=females$mated, n.boot=10000, standardized=TRUE)
f.fem <- moments.differentials(z=females$fem,W=females$mated, n.boot=10000, standardized=TRUE)
f.tib <- moments.differentials(z=females$tib,W=females$mated, n.boot=10000, standardized=TRUE)

# plot abdoment length
f.body <- fitness.landscape(mod=female.gam,phenotype="total_body",covariates=c("tib","fem"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.body$points[,1],f.body$Wbar,type="l", ylim=c(0,1),xlab="Body length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.body$points[,1],f.body$WbarPI[1,],lty=2)
lines(f.body$points[,1],f.body$WbarPI[2,],lty=2)

# plot femur length
f.femur <- fitness.landscape(mod=female.gam,phenotype="fem",covariates=c("tib","total_body"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.femur$points[,1],f.femur$Wbar,type="l", ylim=c(0,1),xlab="Femur length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.femur$points[,1],f.femur$WbarPI[1,],lty=2)
lines(f.femur$points[,1],f.femur$WbarPI[2,],lty=2)

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




# =========================
# male-male competition
# =========================

# dyadic effects

male_male <- male_male
  
library(dplyr)

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

combat2 <- combat2 %>%
  mutate(
    body_diff = total_body_win - total_body_lose,
    fem_diff  = fem_win - fem_lose
  ) # positive = winner is bigger

df <- combat2

set.seed(123)

flip <- rbinom(nrow(df), 1, 0.5)

df$id_i <- ifelse(flip == 1, df$winner_id, df$loser_id)
df$id_j <- ifelse(flip == 1, df$loser_id, df$winner_id)

df$fem_i <- ifelse(flip == 1, df$fem_win, df$fem_lose)
df$fem_j <- ifelse(flip == 1, df$fem_lose, df$fem_win)

df$body_i <- ifelse(flip == 1, df$total_body_win, df$total_body_lose)
df$body_j <- ifelse(flip == 1, df$total_body_lose, df$total_body_win)

df$win <- ifelse(df$id_i == df$winner_id, 1, 0)

df$diff_fem  <- df$fem_i - df$fem_j
df$diff_body <- df$body_i - df$body_j

model <- glm(win ~ diff_fem + diff_body,
             family = binomial, data = df)
summary(model)

# selection analysis of contest outcome

winners <- data.frame(
  id   = combat2$winner_id,
  fem  = combat2$fem_win,
  body = combat2$total_body_win,
  win  = 1
)

losers <- data.frame(
  id   = combat2$loser_id,
  fem  = combat2$fem_lose,
  body = combat2$total_body_lose,
  win  = 0
)

indiv <- rbind(winners, losers)

table(indiv$id)

indiv$fem_z  <- scale(indiv$fem)
indiv$body_z <- scale(indiv$body)

model_sel <- glm(win ~ fem_z + body_z,
                 family = binomial, data = indiv)
summary(model_sel)


## ---- end