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

morph_data <- morph_data %>%
  mutate(fem=rowMeans(select(., l_fem, r_fem), na.rm = TRUE)) %>%
  mutate(tib=rowMeans(select(., l_tib, r_tib), na.rm = TRUE))

# =============================
# selection analysis
# =============================

# male selection analysis
males <- morph_data %>%
  dplyr::filter(sex=="m")

male.gam<-gam(mated~s(tot_abdo)+s(fem)+s(tib),family=binomial(link="logit"),method="GCV.Cp", data=males)
male.grad<-gam.gradients(mod = male.gam,phenotype=c("tot_abdo","fem","tib"),se.method="boot.para",n.boot=1000,standardize=TRUE)

male.estimates<-data.frame(male.grad$ests)

## need to double the (quadratic) gamma gradient
female.estimates[6,1]<-female.estimates[6,1]*2
female.estimates[7,1]<-female.estimates[7,1]*2
female.estimates[8,1]<-female.estimates[8,1]*2
female.estimates[9,1]<-female.estimates[9,1]*2
female.estimates[10,1]<-female.estimates[10,1]*2
save(female.estimates, file = "female.estimates.RData")

#male selection differentials
m.abdomen <- moments.differentials(z=males$tot_abdo,W=males$mated, n.boot=10000, standardized=TRUE)
m.fem <- moments.differentials(z=males$fem,W=males$mated, n.boot=10000, standardized=TRUE)
m.tib <- moments.differentials(z=males$tib,W=males$mated, n.boot=10000, standardized=TRUE)

# plot abdoment length
m.abdomen <- fitness.landscape(mod=male.gam,phenotype="tot_abdo",covariates=c("tib","fem"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(m.abdomen$points[,1],m.abdomen$Wbar,type="l", ylim=c(0,1),xlab="Abdomen length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(m.abdomen$points[,1],m.abdomen$WbarPI[1,],lty=2)
lines(m.abdomen$points[,1],m.abdomen$WbarPI[2,],lty=2)

# plot femur length
m.femur <- fitness.landscape(mod=male.gam,phenotype="fem",covariates=c("tib","tot_abdo"),PI.method="boot.para")
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
vis.gam(male.gam,view=c("tot_abdo","fem"),color="heat",plot.type = "persp",  type="response",theta=155,
        xlab="\nAbdomen length", ylab="\nFemur length",zlab="\nFitness")
male.proj<-recordPlot()

# female selection analysis
females <- morph_data %>%
  dplyr::filter(sex=="f")

female.gam<-gam(mated~s(tot_abdo)+s(fem)+s(tib),family=binomial(link="logit"),method="GCV.Cp", data=females)
female.grad<-gam.gradients(mod = female.gam, phenotype=c("tot_abdo","fem","tib"),se.method="boot.para",n.boot=1000,standardize=TRUE)

female.estimates<-data.frame(female.grad$ests)

## need to double the (quadratic) gamma gradient
female.estimates[6,1]<-female.estimates[6,1]*2
female.estimates[7,1]<-female.estimates[7,1]*2
female.estimates[8,1]<-female.estimates[8,1]*2
female.estimates[9,1]<-female.estimates[9,1]*2
female.estimates[10,1]<-female.estimates[10,1]*2
save(female.estimates, file = "female.estimates.RData")

# female selection differentials
f.abdomen <- moments.differentials(z=females$tot_abdo,W=females$mated, n.boot=10000, standardized=TRUE)
f.fem <- moments.differentials(z=females$fem,W=females$mated, n.boot=10000, standardized=TRUE)
f.tib <- moments.differentials(z=females$tib,W=females$mated, n.boot=10000, standardized=TRUE)

# plot abdoment length
f.abdomen <- fitness.landscape(mod=female.gam,phenotype="tot_abdo",covariates=c("tib","fem"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.abdomen$points[,1],f.abdomen$Wbar,type="l", ylim=c(0,1),xlab="Abdomen length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.abdomen$points[,1],f.abdomen$WbarPI[1,],lty=2)
lines(f.abdomen$points[,1],f.abdomen$WbarPI[2,],lty=2)

# plot femur length
f.femur <- fitness.landscape(mod=female.gam,phenotype="fem",covariates=c("tib","tot_abdo"),PI.method="boot.para")
par(mar=c(6,6,4,4))
plot(f.femur$points[,1],f.femur$Wbar,type="l", ylim=c(0,1),xlab="Femur length (mm)",ylab="Mating success",cex.lab=2,cex.axis=1.5)
lines(f.femur$points[,1],f.femur$WbarPI[1,],lty=2)
lines(f.femur$points[,1],f.femur$WbarPI[2,],lty=2)

dev.off()

# =========================
# female fecundity
# =========================

female.fecundity <- glm.nb(eggs ~ log(tot_abdo), data=females)
summary(female.fecundity) # strong fecundity selection

newdat <- data.frame(
  tot_abdo = seq(min(females$tot_abdo),
                 max(females$tot_abdo),
                 length.out = 100)
)

pred <- predict(female.fecundity, newdata = newdat,
                type = "link", se.fit = TRUE)

newdat$fit <- exp(pred$fit)
newdat$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
newdat$upr <- exp(pred$fit + 1.96 * pred$se.fit)

female.fecundity.plot <- ggplot(females, aes(x = tot_abdo, y = eggs)) +
      geom_point(alpha = 0.5) +
      geom_line(data = newdat, aes(x = tot_abdo, y = fit),
            color = "blue", linewidth = 1) +
      geom_ribbon(data = newdat,
              aes(x = tot_abdo, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              alpha = 0.2) +
      scale_x_log10() +
      labs(x = "Abdomen size (log scale)",
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





## ---- end