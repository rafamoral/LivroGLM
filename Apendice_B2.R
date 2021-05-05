library(tidyverse)
library(ggplot2)

rotenone <- tibble(dose = c(0,2.6,3.8,5.1,7.7,10.2),
                   y = c(0,6,16,24,42,44),
                   m = c(49,50,48,46,49,50))

# Modelo nulo
fit0 <- glm(cbind(y, m - y) ~ 1, family = binomial, data = rotenone)
gl.0 <- fit0$df.residual
X2.0 <- sum(resid(fit0, type = "pearson")^2)
dev0 <- deviance(fit0)

# Modelo de regressão linear
fit1 <- glm(cbind(y, m - y) ~ dose, family = binomial, data = rotenone)
gl.1 <- fit1$df.residual
X2.1 <- sum(resid(fit1, type = "pearson")^2)
dev1 <- deviance(fit1)

Tabela4.3 <- data.frame("Modelo" = c("Nulo", "Reg. linear"),
                        "g.l." = c(gl.0, gl.1),
                        "Deviance" = c(dev0, dev1),
                        "X2" = c(X2.0, X2.1))

# Valores da Tabela 4.4  
anova(fit1, test = "Chisq")

# Figura 4.1
fit1.probit <- update(fit1, family = binomial(link = "probit"))
fit1.cloglog <- update(fit1, family = binomial(link = "cloglog"))

rotenone_pred <- tibble(x = seq(0, 10, length = 30),
                        Logit = predict(fit1, data.frame(dose = x), type = "response"),
                        Probit = predict(fit1.probit, data.frame(dose = x), type = "response"),
                        `Complemento log-log` = predict(fit1.cloglog, data.frame(dose = x), type = "response"))

pdf("Rotenonefit.pdf", w = 7, h = 5)
rotenone_pred %>%
  pivot_longer(2:4,
               names_to = "Função de ligação",
               values_to = "pred") %>%
  ggplot(aes(x = x, y = pred)) +
  theme_bw() +
  geom_point(data = rotenone, aes(x = dose, y = y/m),
             alpha = .5) +
  geom_line(aes(lty = `Função de ligação`)) +
  xlab("Dose") +
  ylab("Proporção de insetos mortos")
dev.off()