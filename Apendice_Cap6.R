########################################################################
## 6.1 Cerejeiras
########################################################################
require(hnp)

## Fig 6.1 - cerejeiras
data(trees, package = "datasets")
names(trees) <- c("D","H","V")

pdf("trees_disp.pdf", w = 9, h = 6)
par(mfrow=c(2,3), cex.lab = 1.2, cex.axis = 1.2)
plot(V ~ D, trees)
plot(V ~ H, trees)
plot(D ~ H, trees)
plot(log(V) ~ log(D), trees)
plot(log(V) ~ log(H), trees)
plot(log(D) ~ log(H), trees)
dev.off()

## Fig. 6.2
mod1 <- lm(V ~ H + D, data = trees)
mod2 <- lm(log(V) ~ log(H) + log(D), data = trees)
mod3 <- glm(V ~ log(H) + log(D), family = gaussian(link = log), data=trees)

par(mfrow = c(1,3), cex.lab = 1.2)
plot(trees$V, fitted(mod1), main = "M1", xlab = "", ylab = "", pch = 16, cex = .7)
lines(loess.smooth(trees$V, fitted(mod1)))
abline(0, 1, lty = 2)
plot(fitted(mod1), rstudent(mod1), main = "M1", xlab = "", ylab = "", pch = 16, cex = .7)
lines(loess.smooth(fitted(mod1), rstudent(mod1)))
abline(h = 0, lty = 2)
hnp(mod1, half = FALSE, xlab = "", ylab = "", main = "M1", resid.type = "student", pch = 16, cex = .7)

par(mfrow = c(1,3), cex.lab = 1.2)
plot(log(trees$V), fitted(mod2), main = "M2", xlab = "", ylab = "Valores ajustados", pch = 16, cex = .7)
lines(loess.smooth(log(trees$V), fitted(mod2)))
abline(0, 1, lty = 2)
plot(fitted(mod2), rstudent(mod2), main = "M2", xlab = "", ylab = "Resíduos", pch = 16, cex = .7)
lines(loess.smooth(fitted(mod2), rstudent(mod2)))
abline(h = 0, lty = 2)
hnp(mod2, half = FALSE, xlab = "", ylab = "Resíduos", main = "M2", resid.type = "student", pch = 16, cex = .7)

par(mfrow = c(1,3), cex.lab = 1.2)
plot(trees$V, fitted(mod3), main = "M3", xlab = "Valores observados", ylab = "", pch = 16, cex = .7)
lines(loess.smooth(trees$V, fitted(mod3)))
abline(0, 1, lty = 2)
plot(fitted(mod3), resid(mod3, type = "deviance"), main = "M3", xlab = "Valores ajustados", ylab = "", pch = 16, cex = .7)
lines(loess.smooth(fitted(mod3), resid(mod3, type = "deviance")))
abline(h = 0, lty = 2)
hnp(mod3, half = FALSE, xlab = "Quantis teóricos", ylab = "", main = "M3", resid.type = "student", pch = 16, cex = .7)

## Figs. 6.3 e 6.4
m1.hnp <- hnp(mod1, half = FALSE, plot = FALSE, resid.type = "student")
m2.hnp <- hnp(mod2, half = FALSE, plot = FALSE, resid.type = "student")
m3.hnp <- hnp(mod3, half = FALSE, plot = FALSE, resid.type = "student")

pdf("fig63.pdf", w = 9, h = 9)
par(mfrow = c(2,2), cex.lab = 1.4, cex.axis = 1.2)
plot(trees$V, fitted(mod1), main = "", xlab = "Valores observados de volumes", ylab = "Valores ajustados", pch = 16, cex = .7)
text(trees$V[c(2,31)], fitted(mod1)[c(2,31)], c("1,2,3","31"), pos = c(4,1))
lines(loess.smooth(trees$V, fitted(mod1)))
abline(0, 1, lty = 2)
plot(abs(dffits(mod1)), main = "", xlab = "Índices", ylab = "Valores absolutos de DFFitS", pch = 16, cex = .7, ylim = c(0,1.5))
abline(h = 3*sqrt(length(coef(mod1))/nrow(trees)), lty=2)
text(31, abs(dffits(mod1))[31], "31", pos = 1)
plot(m1.hnp, pch = 16, xlab = "Quantis teóricos", ylab = "Resíduos estudentizados")
text(m1.hnp$x[c(1,31)], m1.hnp$residuals[c(1,31)], names(m1.hnp$residuals)[c(1,31)], pos = 1)
with(trees, boxcox(V ~ H + D, ylab = "log(função de verossimilhança)"))
dev.off()

pdf("fig64.pdf", w = 9, h = 9)
par(mfrow = c(2,2), cex.lab = 1.4, cex.axis = 1.2)
plot(log(trees$V), fitted(mod2), main = "", xlab = "log(Valores observados de volumes)", ylab = "Valores ajustados", pch = 16, cex = .7)
lines(loess.smooth(log(trees$V), fitted(mod2)))
text(log(trees$V)[c(2,31)], fitted(mod2)[c(2,31)], c("1,2,3","31"), pos = c(4,2))
abline(0, 1, lty = 2)
plot(abs(dffits(mod2)), main = "", xlab = "Índices", ylab = "Valores absolutos de DFFitS", pch = 16, cex = .7, ylim = c(0,1.5))
text(18, abs(dffits(mod2))[18], "18", pos = 1)
abline(h=3*sqrt(length(coef(mod2))/nrow(trees)), lty = 2)
plot(m2.hnp, pch = 16, xlab = "Quantis teóricos", ylab = "Resíduos estudentizados")
text(m2.hnp$x[c(1,31)], m2.hnp$residuals[c(1,31)], names(m2.hnp$residuals)[c(1,31)], pos = c(1,3))
with(trees, boxcox(log(V) ~ log(H) + log(D), ylab = "log(função de verossimilhança)"))
dev.off()

pdf("fig65.pdf", w = 9, h = 9)
par(mfrow = c(2,2), cex.lab = 1.4, cex.axis = 1.2)
plot(trees$V, fitted(mod3), main = "", xlab = "log(Valores observados de volumes)", ylab = "Valores ajustados", pch = 16, cex = .7)
lines(loess.smooth(trees$V, fitted(mod3)))
text(trees$V[c(2,31)], fitted(mod3)[c(2,31)], c("1,2,3","31"), pos = c(4,2))
abline(0, 1, lty = 2)
plot(abs(dffits(mod3)), main = "", xlab = "Índices", ylab = "Valores absolutos de DFFitS", pch = 16, cex = .7, ylim = c(0,1.5))
text(18, abs(dffits(mod3))[18], "18", pos = 1)
abline(h=3*sqrt(length(coef(mod3))/nrow(trees)), lty = 2)
plot(m3.hnp, pch = 16, xlab = "Quantis teóricos", ylab = "Resíduos estudentizados")
text(m3.hnp$x[c(1,31)], m3.hnp$residuals[c(1,31)], names(m3.hnp$residuals)[c(1,31)], pos = c(1,3))
with(trees, boxcox(log(V) ~ log(H) + log(D), ylab = "log(função de verossimilhança)"))
dev.off()

########################################################################
## 6.2 Gordura no leite
########################################################################
library(tidyverse)
library(ggplot2)
require(hnp)

leite <- tibble(gordura = c(0.31,0.39,0.50,0.58,0.59,0.64,0.68,0.66,0.67,0.70,0.72,0.68,
                            0.65,0.64,0.57,0.48,0.46,0.45,0.31,0.33,0.36,0.30,0.26,0.34,
                            0.29,0.31,0.29,0.20,0.15,0.18,0.11,0.07,0.06,0.01,0.01),
                semana = 1:35)

## Modelo normal para log(gordura)
mod1 <- lm(log(gordura) ~ semana + log(semana), data = leite)
summary(mod1)
anova(mod1)

plot(leite$semana, leite$gordura, xlab = "Semanas", ylim = c(0,.9), type = "n",
     ylab = "Produção de gordura (kg/dia) no leite")
x.grid <- seq(1, 35, length = 100)
lines(x.grid, exp(predict(mod1, data.frame(semana = x.grid))))
points(leite$semana, leite$gordura, pch = 21, bg = "lightgray", cex = .85, data = leite)

hnp(mod1, pch = 16, xlab = "Quantis teóricos", half = FALSE,
    ylab = "Resíduos estudentizados", resid.type = "student")

## Modelo normal com função de ligação logarítmica
mod2 <- glm(gordura ~ semana + log(semana),
            family = gaussian(link = "log"),
            data = leite)
summary(mod2)
anova(mod2, test = "F")

plot(leite$semana, leite$gordura, xlab = "Semanas", ylim=c(0,.9), type = "n",
     ylab = "Produção de gordura (kg/dia) no leite")
x.grid <- seq(1, 35, length = 100)
lines(x.grid, exp(predict(mod2, data.frame(semana = x.grid))))
points(leite$semana, leite$gordura, pch = 21, bg = "lightgray", cex = .85)

dfun <- function(obj) rstudent(obj)
sfun <- function(n, obj) {
  dp <- sqrt(summary(obj)$deviance/summary(obj)$df.residual) 
  rnorm(n, obj$fit, dp)
}
ffun <- function(y.) glm(y. ~ semana + log(semana),
                         family = gaussian(link = "log"),
                         start = coef(mod2),
                         data = leite)
hnp(mod2, newclass = TRUE, diagfun=dfun, simfun=sfun, fitfun=ffun,
    pch = 16, xlab = "Quantis teóricos", half = FALSE,
    ylab = "Resíduos estudentizados")

## Gráficos
leite_pred <- tibble(x = seq(1, 35, length = 100),
                     `Transformação logarítmica` = exp(predict(mod1, data.frame(semana = x))),
                     `Ligação logarítmica` = predict(mod2, data.frame(semana = x), type = "response"))

pdf("fat.pdf", w = 7, h = 5)
leite_pred %>%
  pivot_longer(2:3,
               names_to = "modelo",
               values_to = "pred") %>%
  ggplot(aes(x = x, y = pred)) +
  theme_bw() +
  geom_point(data = leite, aes(x = semana, y = gordura),
             alpha = .5) +
  geom_line(aes(lty = modelo)) +
  xlab("Semanas") +
  ylab("Produção de gordura (kg/dia) no leite")
dev.off()

pdf("fat_hnp.pdf", w = 12, h = 6)
par(mfrow = c(1,2), cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4)
hnp(mod1, pch = 16, xlab = "Quantis teóricos", half = FALSE,
    ylab = "Resíduos estudentizados", resid.type = "student", main = "(a)")
hnp(mod2, newclass = TRUE, diagfun = dfun, simfun = sfun, fitfun = ffun,
    pch = 16, xlab = "Quantis teóricos", half = FALSE,
    ylab = "Resíduos estudentizados", main = "(b)")
dev.off()

########################################################################
## 6.3 Acácia Negra
########################################################################

require(hnp)
require(tidyverse)
require(gridExtra)

acacia <- read.csv2("dados_volume.csv", header = TRUE)

# análise exploratória
ac1 <- acacia %>%
  ggplot(aes(x = d, y = v)) +
  theme_bw() +
  facet_wrap(~ local) +
  geom_point() +
  xlab("DAP") +
  ylab("Volume de madeira")

ac2 <- acacia %>%
  ggplot(aes(x = h, y = v)) +
  theme_bw() +
  facet_wrap(~ local) +
  geom_point() +
  xlab("Altura") +
  ylab("Volume de madeira")

ac3 <- acacia %>%
  ggplot(aes(x = h, y = d)) +
  theme_bw() +
  facet_wrap(~ local) +
  geom_point() +
  xlab("Altura") +
  ylab("DAP")

pdf("acacia_exp.pdf", w = 10, h = 10)
grid.arrange(ac1, ac2, ac3)
dev.off()

# ajuste do modelo normal com função de ligação logarítmica
mod_normal <- glm(v ~ local * (log(h) + log(d)),
                  family = gaussian(link = log),
                  data = acacia)

dfun <- function(obj) rstudent(obj)
sfun <- function(n, obj) simulate(obj)$sim_1
ffun <- function(resp) glm(resp ~ local * (log(h) + log(d)),
                           family = gaussian(link = log),
                           data = acacia,
                           start = coef(mod_normal))

set.seed(2022)
hnp_normal <- hnp(mod_normal, newclass = TRUE,
                  diagfun = dfun, simfun = sfun, fitfun = ffun)

# ajuste do modelo normal com função de ligação logarítmica
mod_gama <- glm(v ~ local * (log(h) + log(d)),
                family = Gamma(link = "log"), 
                data = acacia)

set.seed(2022)
hnp_gama <- hnp(mod_gama, resid.type = "deviance")

pdf("acacia_hnp.pdf", w = 12, h = 6)
par(mfrow = c(1,2))
plot(hnp_normal, pch = 16, xlab = "Quantis da distribuição meio-normal",
     ylab = "Resíduos estudentizados", main = "(a)")
plot(hnp_gama, pch = 16, xlab = "Quantis da distribuição meio-normal",
     ylab = "Resíduos componentes do desvio", main = "(b)")
dev.off()

# testando interações duplas
drop1(mod_gama, test = "F")

# testando efeitos principais
mod_gama2 <- glm(v ~ local + log(h) + log(d),
                 family = Gamma(link = "log"), 
                 data = acacia)
drop1(mod_gama2, test = "F")

# gráficos de diagnósticos
pdf("acacia_diag.pdf", w = 9, h = 9)
par(mfrow = c(2,2))
plot(acacia$v, fitted(mod_normal), main = "Modelo normal", xlab = "Valores observados", ylab = "Valores ajustados", pch = 16, cex = .7)
lines(loess.smooth(acacia$v, fitted(mod_normal)))
abline(0, 1, lty = 2)
plot(fitted(mod_normal), rstudent(mod_normal), main = "Modelo normal", xlab = "Valores ajustados", ylab = "Resíduos estudentizados", pch = 16, cex = .7)
lines(loess.smooth(fitted(mod_normal), rstudent(mod_normal)))
abline(h = 0, lty = 2)
plot(acacia$v, fitted(mod_gama), main = "Modelo gama", xlab = "Valores observados", ylab = "Valores ajustados", pch = 16, cex = .7)
lines(loess.smooth(acacia$v, fitted(mod_gama)))
abline(0, 1, lty = 2)
plot(fitted(mod_gama), resid(mod_gama, type = "deviance"), main = "Modelo gama", xlab = "Valores ajustados", ylab = "Resíduos componentes do desvio", pch = 16, cex = .7)
lines(loess.smooth(fitted(mod_gama), resid(mod_gama, type = "deviance")))
abline(h = 0, lty = 2)
dev.off()

########################################################################
## 6.4 Sobrevivência de ratos
########################################################################
library(tidyverse)
library(ggplot2)
require(hnp)

rato.dat <- read.table("rato.txt", header = TRUE) %>% as_tibble
rato.dat$tipo <- as.factor(rato.dat$tipo)
rato.dat$trat <- as.factor(rato.dat$trat)

plot(y ~ tipo : trat, rato.dat)

rato.dat <- rato.dat %>%
  mutate(y_inv = 1/y,
         y_3_4 = y^(-3/4))

pdf("Fig-Rato-boxplot.pdf", w = 12, h = 4)
rato.dat %>%
  pivot_longer(c(1,4,5),
               names_to = "transformacao",
               values_to = "y") %>%
  mutate(transformacao = recode_factor(transformacao, "y" = "Y", "y_inv" = "1/Y", "y_3_4" = "Y^(-3/4)")) %>%
  ggplot(aes(y = y, x = tipo : trat)) +
  theme_bw() +
  geom_boxplot() +
  facet_wrap(~ transformacao) +
  xlab("Tipos de venenos e tratamentos") +
  ylab("Tempo de sobrevivência")
dev.off()

pdf("Fig-Rato-boxcox.pdf", w = 12, h = 4)
par(mfrow = c(1,3), cex.lab = 1.5, cex.axis = 1.3)
boxcox(y ~ tipo * trat, data = rato.dat, 
       ylab = "Log(função de verossimilhança)")
title(expression(Y), cex.main = 2.2)
boxcox(1/y ~ tipo * trat, data = rato.dat, 
       ylab = "Log(função de verossimilhança)")
title(expression(1/Y), cex.main = 2.2)
boxcox(y^{-3/4} ~ tipo * trat, data = rato.dat, 
       ylab = "Log(função de verossimilhança)")
title(expression(Y^{-3/4}), cex.main = 2.2)
dev.off()

## Ajuste dos modelos
mod1 <- lm(y ~ tipo * trat, data = rato.dat)
mod2 <- lm(1/y ~ tipo * trat, data = rato.dat)
mod3 <- lm(y^{-3/4} ~ tipo * trat, data = rato.dat)

m1hnp <- hnp(mod1, half = FALSE)
m2hnp <- hnp(mod2, half = FALSE)
m3hnp <- hnp(mod3, half = FALSE)

pdf("Fig-Rato-todos.pdf", w = 12, h = 12)
par(mfrow = c(3,3), cex = .8, cex.main = 2, cex.lab = 1.5, cex.axis = 1.3)
plot(fitted(mod1) ~ rato.dat$y, xlab = "Valores observados",
     ylab = "Valores ajustados", main = expression(Y))
plot(fitted(mod2) ~ 1/rato.dat$y, xlab = "1/(Valores observados)",
     ylab = "Valores ajustados", main = expression(1/Y))
plot(fitted(mod3) ~ I(rato.dat$y^{-3/4}), xlab = "(Valores observados)^(-3/4)",
     ylab = "Valores ajustados", main = expression(Y^{-3/4}))
plot(rstudent(mod1) ~ fitted(mod1), xlab = "Valores ajustados",
     ylab = "Resíduos"); abline(h = 0, lty = 2)
plot(rstudent(mod2) ~ fitted(mod2), xlab = "Valores ajustados",
     ylab = "Resíduos"); abline(h = 0, lty = 2)
plot(rstudent(mod3) ~ fitted(mod3), xlab = "Valores ajustados",
     ylab = "Resíduos"); abline(h = 0, lty = 2)
plot(m1hnp, pch = 16, xlab = "Resíduos", ylab = "Quantis teóricos", cex = .85)
plot(m2hnp, pch = 16, xlab = "Resíduos", ylab = "Quantis teóricos", cex = .85)
plot(m3hnp, pch = 16, xlab = "Resíduos", ylab = "Quantis teóricos", cex = .85)
dev.off()

########################################################################
## 6.5 Assinaturas de TV a cabo
########################################################################

require(hnp)

tvcabo <- read.table("tv-cabo.txt", header = TRUE)
# y = número de assinantes de TV a cabo (em milhares)
# x1 = número de domicílios na área
# x2 = renda per capita por domicílio
# x3 = taxa de instalação
# x4 = custo médio mensal de manutenção
# x5 = número de canais disponíveis
# x6 = número de canais não pagos

# Gráficos de dispersão
pdf("Fig-tv-disp.pdf", w = 9, h = 12)
par(mfrow = c(4,3), cex.main = 2, cex.lab = 1.6, cex.axis = 1.4)
plot(y ~ x1, data = tvcabo, xlab = expression(x[1]))
plot(y ~ x2, data = tvcabo, xlab = expression(x[2]))
plot(y ~ x3, data = tvcabo, xlab = expression(x[3]))
plot(y ~ x4, data = tvcabo, xlab = expression(x[4]))
plot(y ~ x5, data = tvcabo, xlab = expression(x[5]))
plot(y ~ x6, data = tvcabo, xlab = expression(x[6]))
plot(y ~ log(x1), data = tvcabo, xlab = expression(log(x[1])))
plot(y ~ log(x2), data = tvcabo, xlab = expression(log(x[2])))
plot(y ~ log(x3), data = tvcabo, xlab = expression(log(x[3])))
plot(y ~ log(x4), data = tvcabo, xlab = expression(log(x[4])))
plot(y ~ log(x5), data = tvcabo, xlab = expression(log(x[5])))
plot(y ~ log(x6), data = tvcabo, xlab = expression(log(x[6])))
dev.off()

# Potência máxima de Box-Cox
pdf("Fig-tv-boxcox.pdf", w = 6, h = 6)
par(mfrow = c(1,1), cex.lab = 1.6, cex.axis = 1.4)
with(tvcabo, boxcox(y ~ log(x1) + log(x2) + log(x3) + log(x4) + log(x5) + log(x6),
                    ylab = "Log(função de verossimilhança"))
dev.off()

# Modelo M1
m1 <- lm(y ~ log(x1) + log(x2) + log(x3) + log(x4) + log(x5) + log(x6),
         data = tvcabo)
summary(m1)

# Modelo M2
m2 <- lm(log(y) ~ log(x1) + log(x2) + log(x3) + log(x4) + log(x5) + log(x6),
         data = tvcabo)
summary(m2)
dropterm(m2, test = "F")

m2semx3 <- lm(log(y) ~ log(x1) + log(x2) + log(x4) + log(x5) + log(x6),
             data = tvcabo)
summary(m2semx3)

## Gráficos
m2semx3hnp <- hnp(m2semx3, half = F)
xx <- m2semx3hnp$x
names(xx) <- names(m2semx3hnp$residuals)

pdf("Fig-tv-fit.pdf", w = 9, h = 3)
par(mfrow = c(1,3), cex.lab = 1.5, cex.axis = 1.3, cex.main = 2)
plot(fitted(m2semx3) ~ log(tvcabo$y), xlab = "Log(Valores observados)",
     ylab = "Valores ajustados", main = "(a)", cex = .9)
text(log(tvcabo$y)[c(11,14,26)], fitted(m2semx3)[c(11,14,26)],
     c(11,14,26), cex = 1.2, pos = c(1,1,3))
abline(0, 1, lty = 2)
plot(rstudent(m2semx3) ~ fitted(m2semx3), xlab = "Valores ajustados",
     ylab = "Resíduos", main = "(b)", cex = .9)
text(fitted(m2semx3)[c(11,14,26)], rstudent(m2semx3)[c(11,14,26)],
     c(11,14,26), cex = 1.2, pos = c(3,1,3))
abline(h = 0, lty = 2)
plot(m2semx3hnp, xlab = "Quantis teóricos", ylab = "Resíduos",
     pch = 16, main = "(c)", cex = .9)
text(xx[c("11","14","26")], m2semx3hnp$residuals[c("11","14","26")],
     c(11,14,26), cex = 1.2, pos = c(4,3,1))
dev.off()