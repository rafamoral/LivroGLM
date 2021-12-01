########################################################################
## 7.1 Rotenone
########################################################################
library(tidyverse)
library(ggplot2)

rotenone <- data.frame(dose = c(0,2.6,3.8,5.1,7.7,10.2),
                       y = c(0,6,16,24,42,44),
                       m = c(49,50,48,46,49,50))

fit <- glm(cbind(y, m-y) ~ dose, family = binomial, data = rotenone)
# Dose letal 50
(d50 <- - coef(fit)[1]/coef(fit)[2])

# Encontrando a dose letal utilizando o pacote MASS
require(MASS)
(v <- vcov(fit)) # matriz de variâncias e covariâncias

dose.p(fit) ## p=0.5
dose.p(fit, p = 0.9)  ## p=0.9

# Doses LD25, LD50, LD75
dose.p(fit, p = 1:3/4)

# Obtendo os intervalos de confiança para a DL50
source("ic_dose_letal.R")
confint(dose.p(fit))
Fieller(fit)
# gráfico do perfil de deviance
LR_confint(fit, profile = TRUE, xlab = expression(theta))

########################################################################
## 7.2 Cipermetrina
########################################################################
# Collett (1991), página 75
## Pacotes necessários
require(hnp)
require(latticeExtra)

## Função para extrair a estatística X2
X2 <- function(obj) sum(resid(obj, type="pearson")^2)

## Carregando conjunto de dados
cyper <- data.frame(y = c(1,4,9,13,18,20,0,2,6,10,12,16),
                    m = 20,
                    dose = c(1,2,4,8,16,32),
                    sexo = gl(2, 6, labels=c("Macho","Fêmea")))
cyper$ldose <- log(cyper$dose, 2)

## Gráfico exploratório
cyper %>%
  ggplot(aes(x = log(dose, 2), y = y/m)) +
  theme_bw() +
  geom_point(aes(pch = sexo)) +
  geom_line(aes(lty = sexo)) +
  xlab("log(dose)") +
  ylab("Proporção de insetos mortos")

## Ajustando modelos utilizando dose como fator qualitativo
mod1 <- glm(cbind(y, m-y) ~ 1, family = binomial, data = cyper) # constante
mod2 <- glm(cbind(y, m-y) ~ sexo, family = binomial, data = cyper) # sexo
mod3 <- glm(cbind(y, m-y) ~ factor(dose), family = binomial, data = cyper) # dose
mod4 <- glm(cbind(y, m-y) ~ sexo + factor(dose), family = binomial, data = cyper) # sexo + dose|sexo
mod5 <- glm(cbind(y, m-y) ~ factor(dose) + sexo, family = binomial, data = cyper) # dose + sexo|dose

## Tabela 7.2
Modelo <- g.l. <- Desvios <- X.2 <- NULL
for(i in 1:4) {
  obj <- get(paste("mod", i, sep = ""))
  Modelo[i] <- as.character(formula(obj))[3]
  g.l.[i] <- obj$df.residual
  Desvios[i] <- deviance(obj)
  X.2[i] <- X2(obj)
}

tabela7.2 <- data.frame(Modelo, g.l., Desvios = round(Desvios, 2), 
                        p = round(pchisq(Desvios, g.l., lower = FALSE), 4), X.2 = round(X.2, 2),
                        p. = round(pchisq(X.2, g.l., lower = FALSE), 4))
tabela7.2

## Tabela 7.3
s.d <- anova(mod4, test = "Chisq")
d.s <- anova(mod5, test = "Chisq")
tabela7.3 <- rbind(rbind(s.d, d.s)[c(2,5,6,3),-(3:4)],
                   c("Df" = s.d[3,3], "Deviance" = s.d[3,4],
                     "Pr(>Chi)" = pchisq(s.d[3,4], s.d[3,3], lower = FALSE)),
                   c("Df" = mod4$df.null, "Deviance" = mod4$null.deviance, NA))
tabela7.3$Deviance <- round(tabela7.3$Deviance, 2)
tabela7.3[,3] <- round(tabela7.3[,3], 4)
row.names(tabela7.3) <- c("Sexo","Sexo|Dose","Dose","Dose|Sexo","Resíduo","Total")
tabela7.3

## Ajustando modelos utilizando dose como fator quantitativo (na escala log)
m1 <- glm(cbind(y, m-y) ~ 1, family = binomial, data = cyper) # modelo nulo
m2 <- glm(cbind(y, m-y) ~ sexo, family = binomial, data = cyper) # retas paralelas ao eixo x
m3 <- glm(cbind(y, m-y) ~ ldose, family = binomial, data = cyper) # retas coincidentes
m4 <- glm(cbind(y, m-y) ~ sexo + ldose, family = binomial, data = cyper) # retas paralelas
m5 <- glm(cbind(y, m-y) ~ ldose / sexo, family = binomial, data = cyper) # retas com intercepto comum
m6 <- glm(cbind(y, m-y) ~ ldose * sexo, family = binomial, data = cyper) # retas concorrentes

## Testes para modelos encaixados
anova(m1, m2, m4, m6, test = "Chisq")
anova(m1, m3, m4, m6, test = "Chisq")
anova(m1, m3, m5, m6, test = "Chisq")

## Tabela 7.4
Modelo <- g.l. <- Desvios <- X.2 <- NULL
for(i in c(1,3,4,5,6)) {
  obj <- get(paste("m", i, sep = ""))
  Modelo[i] <- as.character(formula(obj))[3]
  g.l.[i] <- obj$df.residual
  Desvios[i] <- deviance(obj)
  X.2[i] <- X2(obj)
}

tabela7.4 <- data.frame(Modelo, g.l., Desvios = round(Desvios, 2), 
                        p = round(pchisq(Desvios, g.l., lower = FALSE), 4), X.2 = round(X.2, 2),
                        p. = round(pchisq(X.2, g.l., lower = FALSE), 4))[-2,]
tabela7.4

## Tabela 7.5
am4 <- anova(m4, test = "Chisq")
tabela7.5 <- rbind(am4[c(2,3),-(3:4)],
                   c("Df" = am4[3,3], "Deviance" = am4[3,4],
                     "Pr(>Chi)" = pchisq(am4[3,4], am4[3,3], lower = FALSE)),
                   c("Df" = mod4$df.null, "Deviance" = mod4$null.deviance, NA))
tabela7.5$Deviance <- round(tabela7.5$Deviance, 2)
tabela7.5[,3] <- round(tabela7.5[,3], 4)
row.names(tabela7.5) <- c("Sexo","Regressão linear","Resíduo","Total")
tabela7.5

## Doses letais (LD50)
require(MASS)
m4.2 <- glm(cbind(y, m-y) ~ sexo + ldose - 1, family = binomial, data = cyper)
v <- vcov(m4.2)
coefi <- coef(m4.2)

## Utilizando o pacote rootSolve
require(rootSolve)
uniroot.all(males <- function(x, p = 0.5) {
  coefi[1] + coefi[3]*log(x, 2) - log(p/(1-p))
}, c(0, 10))

uniroot.all(females <- function(x, p = 0.5) {
  coefi[2] + coefi[3]*log(x, 2) - log(p/(1-p))
}, c(0, 10))

## Utilizando -beta0/beta1
thetam <- -coefi[1]/coefi[3]
thetaf <- -coefi[2]/coefi[3]
2^thetaf
2^thetam

## Erros-padrões
sqrt(sigma2m <- (v[1,1] + thetam^2*v[3,3] + 2*thetam*v[1,3])/coefi[3]^2) ## fêmea
sqrt(sigma2f <- (v[2,2] + thetaf^2*v[3,3] + 2*thetaf*v[2,3])/coefi[3]^2) ## macho

source("ic_dose_letal.R")
## Método delta
confint(dose.p(m4.2, cf = c(1,3)))
confint(dose.p(m4.2, cf = c(2,3)))
2^confint(dose.p(m4.2, cf = c(1,3)))
2^confint(dose.p(m4.2, cf = c(2,3)))

## Fieller
Fieller(m4.2, cf = c(1,3))
Fieller(m4.2, cf = c(2,3))
2^Fieller(m4.2, cf = c(1,3))
2^Fieller(m4.2, cf = c(2,3))

## Perfil de verossimilhanças
LR_confint(m4.2, cf = c(1,3))
LR_confint(m4.2, cf = c(2,3))
2^LR_confint(m4.2, cf = c(1,3))
2^LR_confint(m4.2, cf = c(2,3))

## Figura 7.1
cyper.pred <- expand.grid(ldose = seq(0, 5, length = 100),
                          sexo = levels(cyper$sexo))
cyper.pred$pred <- predict(m4, cyper.pred, type = "response")

par(mfrow = c(1,2))
hnp(m4, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", main = "(a)", pch = 16, cex = .7)
plot(y/m ~  log(dose, 2), data = cyper, col = 1, main = "(b)",
     pch = rep(16:17, each = 6), xlab = expression(log[2](dose)), 
     ylab = "Proporção de insetos mortos", cex = .8)
lines(pred ~ ldose, lwd=2, data=subset(cyper.pred, cyper.pred$sexo=="Macho"))
lines(pred ~ ldose, lwd=2, lty=2, data=subset(cyper.pred, cyper.pred$sexo=="Fêmea"))
legend("bottomright", c("Machos","Fêmeas"), cex=.8,
       pch=16:17, lty=1:2, col=1, inset=.01, bty="n", lwd=1, y.intersp=1.5)

########################################################################
## 7.3 Mortalidade do besouro da farinha
########################################################################
# Collett(2002) - pág. 103
# Entrada dos dados
tribolium <- data.frame(y = c(3,5,19,19,24,35,2,14,20,27,41,40,28,37,46,48,48,50),
                        m = c(50,49,47,50,49,50,50,49,50,50,50,50,50,50,50,50,50,50),
                        inseticida = gl(3, 6, labels = c("DDT","BHC","mistura")),
                        dose = c(2.00,2.64,3.48,4.59,6.06,8.00))
tribolium$ldose <- log(tribolium$dose)
tribolium$prop <- with(tribolium, y/m)

# Gráficos exploratórios
tribolium %>%
  ggplot(aes(x = dose, y = prop)) +
  theme_bw() +
  geom_point(aes(colour = inseticida)) +
  xlab("Dose") +
  ylab("Proporções de insetos mortos") +
  ggtitle("Usando dose")

tribolium %>%
  ggplot(aes(x = ldose, y = prop)) +
  theme_bw() +
  geom_point(aes(colour = inseticida)) +
  xlab("log(dose)") +
  ylab("Proporções de insetos mortos") +
  ggtitle("Usando log(dose)")

# Ajuste dos modelos usando dose como fator qualitativo (modelo1 -- modelo4)

modelo1 <- glm(cbind(y, m-y) ~ 1, family = binomial, data = tribolium)
modelo2 <- glm(cbind(y, m-y) ~ factor(dose), family = binomial, data = tribolium)
modelo3 <- glm(cbind(y, m-y) ~ inseticida, family = binomial, data = tribolium)
modelo4 <- glm(cbind(y, m-y) ~ factor(dose) + inseticida, family = binomial, data = tribolium)

deviance(modelo1); pchisq(deviance(modelo1), df.residual(modelo1), lower.tail = FALSE)
deviance(modelo2); pchisq(deviance(modelo2), df.residual(modelo2), lower.tail = FALSE)
deviance(modelo3); pchisq(deviance(modelo3), df.residual(modelo3), lower.tail = FALSE)
deviance(modelo4); pchisq(deviance(modelo4), df.residual(modelo4), lower.tail = FALSE)

X2 <- function(obj) sum(resid(obj, type="pearson")^2)
X2(modelo1); pchisq(X2(modelo1), df.residual(modelo1), lower.tail = FALSE)
X2(modelo2); pchisq(X2(modelo2), df.residual(modelo2), lower.tail = FALSE)
X2(modelo3); pchisq(X2(modelo3), df.residual(modelo3), lower.tail = FALSE)
X2(modelo4); pchisq(X2(modelo4), df.residual(modelo4), lower.tail = FALSE)

anova(modelo1, modelo2, modelo4, test="Chisq")
anova(modelo1, modelo3, modelo4, test="Chisq")

# Ajuste dos modelos usando dose como fator quantitativo (mod1 -- mod6)

## Modelo nulo
mod1 <- glm(cbind(y, m-y) ~ 1, family = binomial, data = tribolium)
summary(mod1)
pchisq(deviance(mod1), df.residual(mod1), lower.tail = FALSE)
X2.m1 <- sum(residuals(mod1, "pearson")^2)
pchisq(X2.m1, df.residual(mod1), lower.tail = FALSE)

## Modelo de retas coincidentes
mod2 <- glm(cbind(y, m-y) ~ dose, family = binomial, data=tribolium)
summary(mod2)
pchisq(deviance(mod2), df.residual(mod2), lower.tail = FALSE)
X2.m2 <- sum(residuals(mod2, "pearson")^2)
pchisq(X2.m2, df.residual(mod2), lower.tail = FALSE)

## Modelo de retas paralelas ao eixo x
mod3 <- glm(cbind(y, m-y) ~ inseticida, family = binomial, data = tribolium)
summary(mod3)
pchisq(deviance(mod3), df.residual(mod3), lower.tail = FALSE)
X2.m3 <- sum(residuals(mod3, "pearson")^2)
pchisq(X2.m3, df.residual(mod3), lower.tail = FALSE)

## Modelo de retas paralelas
mod4 <- glm(cbind(y, m-y) ~ inseticida + dose, family = binomial, data = tribolium)
summary(mod4)
pchisq(deviance(mod4), df.residual(mod4), lower.tail = FALSE)
X2.m4 <- sum(residuals(mod4, "pearson")^2)
pchisq(X2.m4, df.residual(mod4), lower.tail = FALSE)
anova(mod4)

## Modelo de retas com intercepto comum
mod5 <- glm(cbind(y, m-y) ~ inseticida : dose, family = binomial, data = tribolium)
summary(mod5)
pchisq(deviance(mod5), df.residual(mod5), lower.tail = FALSE)
X2.m5 <- sum(residuals(mod5, 'pearson')^2)
pchisq(X2.m5, df.residual(mod5), lower.tail = FALSE)

## Modelo de retas concorrentes
mod6 <- glm(cbind(y, m-y) ~ inseticida * dose, family = binomial, data = tribolium)
summary(mod6)
pchisq(deviance(mod6), df.residual(mod6), lower.tail = FALSE)
X2.m6 <- sum(residuals(mod6, "pearson")^2)
pchisq(X2.m6, df.residual(mod6), lower.tail = FALSE)

## Testes para modelos encaixados
anova(mod1, mod2, mod4, mod6, test = "Chisq")
anova(mod1, mod2, mod5, mod6, test = "Chisq")
anova(mod1, mod3, mod4, mod6, test = "Chisq")

# Ajuste dos modelos usando log(dose) (mod7 -- mod10)

## Modelo de retas coincidentes
mod7 <- glm(cbind(y, m-y) ~ ldose, family = binomial, data = tribolium)
summary(mod7)
deviance(mod7)
pchisq(deviance(mod7), df.residual(mod7), lower.tail = FALSE)
X2(mod7)
pchisq(X2(mod7), df.residual(mod7), lower.tail = FALSE)

## Modelo de retas paralelas
mod8 <- glm(cbind(y, m-y) ~ inseticida + ldose, family = binomial, data = tribolium)
summary(mod8)
deviance(mod8)
pchisq(deviance(mod8), df.residual(mod8), lower.tail = FALSE)
X2(mod8)
pchisq(X2(mod8), df.residual(mod8), lower.tail = FALSE)
anova(mod8, test="Chisq")

## Modelo de retas com intercepto comum
mod9 <- glm(cbind(y, m-y) ~ inseticida : ldose, family = binomial, data = tribolium)
summary(mod9)
deviance(mod9)
pchisq(deviance(mod9), df.residual(mod9), lower.tail = FALSE)
X2(mod9)
pchisq(X2(mod9), df.residual(mod9), lower.tail = FALSE)

## Modelo de retas concorrentes
mod10 <- glm(cbind(y, m-y) ~ inseticida * ldose, family = binomial, data = tribolium)
summary(mod10)
deviance(mod10)
pchisq(deviance(mod10), df.residual(mod10), lower.tail = FALSE)
X2(mod10)
pchisq(X2(mod10), df.residual(mod10), lower.tail = FALSE)

## Testes para modelos encaixados
anova(mod1, mod3, mod8, mod10, test = "Chisq")
anova(mod1, mod7, mod8, mod10, test = "Chisq")
anova(mod1, mod7, mod9, mod10, test = "Chisq")

# Half-normal plots com envelope de simulação
require(hnp)
hnp(mod1, print = TRUE)
hnp(mod2, print = TRUE)
hnp(mod3, print = TRUE)
hnp(mod4, print = TRUE)
hnp(mod5, print = TRUE)
hnp(mod6, print = TRUE)
hnp(mod7, print = TRUE)
hnp(mod8, print = TRUE)
hnp(mod9, print = TRUE)
hnp(mod10, print = TRUE)

## Valores preditos

pred <- expand.grid(dose = seq(2, 8, length = 30),
                    inseticida = levels(tribolium$inseticida))
predlog <- expand.grid(ldose = log(seq(2, 8, length = 30)),
                       inseticida = levels(tribolium$inseticida))

pred$pi <- predict(mod4, pred, type = "response")
predlog$pi <- predict(mod8, predlog, type = "response")

## Curvas
par(mfrow = c(1,2))
hnp(mod8, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", main = "(a)", pch = 16, cex = .7)
plot(prop ~ ldose, data = tribolium, col = 1, cex = .8, main = "(b)",
     xlab = "log(dose)", ylab = "Proporções de insetos mortos",
     pch = rep(c(16,1,17), each = 6))
lines(pi ~ ldose, data = subset(predlog, predlog$inseticida == "DDT"), lwd = 2)
lines(pi ~ ldose, data = subset(predlog, predlog$inseticida == "BHC"), lty = 3, lwd = 2)
lines(pi ~ ldose, data = subset(predlog, predlog$inseticida == "mistura"), lty = 4, lwd = 2)
legend("bottomright", c("DDT","BHC","Mistura"), col = 1, pch = c(16,1,17),
       cex = .8, inset = .01, bty = "n", lty = c(1,3,4), y.intersp = 1.5)

# Agrupando DDT e BHC
levels(tribolium$inseticida)
tribolium$inseticida2 <- tribolium$inseticida
levels(tribolium$inseticida2) <- c(1, 1, 2)

mod11 <- glm(cbind(y, m-y) ~ inseticida2 + ldose, family = binomial, data = tribolium)
anova(mod11, mod8, test = "Chisq")

##### Doses letais e intervalos de confiança

## Doses letais
mod8b <- glm(cbind(y, m-y) ~ inseticida + ldose - 1, family = binomial, data = tribolium)
coefi <- coef(mod8b)
dl_DDT <- -coefi[1]/coefi[4]
dl_BHC <- -coefi[2]/coefi[4]
dl_mistura <- -coefi[3]/coefi[4]
exp(dl_DDT)
exp(dl_BHC)
exp(dl_mistura)

source("ic_dose_letal.R")
## Método delta
mod8b <- glm(cbind(y, m-y) ~ inseticida + ldose - 1, family = binomial, data = tribolium)
coef(mod8b)
confint(dose.p(mod8b, cf = c(1,4)))
confint(dose.p(mod8b, cf = c(2,4)))
confint(dose.p(mod8b, cf = c(3,4)))
exp(confint(dose.p(mod8b, cf = c(1,4))))
exp(confint(dose.p(mod8b, cf = c(2,4))))
exp(confint(dose.p(mod8b, cf = c(3,4))))

## Método de Fieller
Fieller(mod8b, cf = c(1,4))
Fieller(mod8b, cf = c(2,4))
Fieller(mod8b, cf = c(3,4))
exp(Fieller(mod8b, cf = c(1,4)))
exp(Fieller(mod8b, cf = c(2,4)))
exp(Fieller(mod8b, cf = c(3,4)))

## Perfil de verossimilhanças
par(mfrow=c(1,3))
LR_confint(mod8b, cf = c(1,4), profile = TRUE)
LR_confint(mod8b, cf = c(2,4), profile = TRUE)
LR_confint(mod8b, cf = c(3,4), profile = TRUE)
exp(LR_confint(mod8b, cf = c(1,4)))
exp(LR_confint(mod8b, cf = c(2,4)))
exp(LR_confint(mod8b, cf = c(3,4)))

########################################################################
## 7.4 Proporções de gemas florais de macieiras
########################################################################
## Pacotes necessários
require(hnp)
X2 <- function(obj) sum(resid(obj, type="pearson")^2)

## Entrada dos dados
gema <- data.frame(variedade = gl(3, 5, labels = c("Crispin","Cox","Golden Delicious")),
                   frutos = 0:4,
                   total.gemas = c(69,93,147,149,151,34,92,133,146,111,21,89,118,124,81),
                   gemas.florais = c(42,43,59,57,43,12,15,18,14,9,6,20,20,21,4))
gema$proporcao <- with(gema, gemas.florais/total.gemas)

resp <- with(gema, cbind(gemas.florais, total.gemas-gemas.florais))

## Gráficos exploratórios
gema %>%
  ggplot(aes(x = frutos, y = proporcao)) +
  theme_bw() +
  geom_point(aes(pch = variedade)) +
  geom_line(aes(lty = variedade)) +
  xlab("Número de frutos") +
  ylab("Proporção de gemas florais")

## Ajustes dos modelos
m1 <- glm(resp ~ 1, family = binomial, data = gema)
deviance(m1)
pchisq(deviance(m1), df.residual(m1), lower.tail = FALSE)
X2(m1)
pchisq(X2(m1), df.residual(m1), lower.tail = FALSE)

m2 <- glm(resp ~ frutos, family = binomial, data = gema)
deviance(m2)
pchisq(deviance(m2), df.residual(m2), lower.tail = FALSE)
X2(m2)
pchisq(X2(m2), df.residual(m2), lower.tail = FALSE)

m3 <- glm(resp ~ variedade, family = binomial, data = gema)
deviance(m3)
pchisq(deviance(m3), df.residual(m3), lower.tail = FALSE)
X2(m3)
pchisq(X2(m3), df.residual(m3), lower.tail = FALSE)

m4 <- glm(resp ~ frutos / variedade, family = binomial, data = gema)
deviance(m4)
pchisq(deviance(m4), df.residual(m4), lower.tail = FALSE)
X2(m4)
pchisq(X2(m4), df.residual(m4), lower.tail = FALSE)

m5 <- glm(resp ~ frutos + variedade, family = binomial, data = gema)
deviance(m5)
pchisq(deviance(m5), df.residual(m5), lower.tail = FALSE)
X2(m5)
pchisq(X2(m5), df.residual(m5), lower.tail = FALSE)

m6 <- glm(resp ~ frutos * variedade, family = binomial, data = gema)
deviance(m6)
pchisq(deviance(m6), df.residual(m6), lower.tail = FALSE)
X2(m6)
pchisq(X2(m6), df.residual(m6), lower.tail = FALSE)

## Testes para modelos encaixados
anova(m1, m2, m4, m6, test = "Chisq")
anova(m1, m2, m5, m6, test = "Chisq")
anova(m1, m3, m5, m6, test = "Chisq")

## Half-normal plots
hnp(m1, print = TRUE)
hnp(m2, print = TRUE)
hnp(m3, print = TRUE)
hnp(m4, print = TRUE)
hnp(m5, print = TRUE)
hnp(m6, print = TRUE)

## Curvas preditas
gema.pred <- expand.grid(frutos = seq(0, 4, length = 100),
                         variedade = levels(gema$variedade))
gema.pred$pred <- predict(m5, gema.pred, type = "response")

gema %>%
  ggplot(aes(x = frutos, y = proporcao)) +
  theme_bw() +
  geom_point(aes(pch = variedade)) +
  geom_line(data = gema.pred, aes(y = pred, lty = variedade)) +
  xlab("Número de frutos") +
  ylab("Proporção de gemas florais")

## Testando diferença entre Cox e Golden Delicious
gema$variedade2 <- gema$variedade
levels(gema$variedade2) <- c(1,2,2)
m5.2 <- glm(resp ~ variedade2 + frutos, family = binomial, data = gema)
anova(m5.2, m5, test = "Chisq")
anova(m5.2, test = "Chisq")
hnp(m5.2, print = TRUE)
coef(update(m5.2, . ~ . - 1))

## Novas curvas preditas
gema.pred2 <- expand.grid(frutos = seq(0, 4, length = 100),
                          variedade2 = levels(gema$variedade2))
gema.pred2$pred <- predict(m5.2, gema.pred2, type = "response")

gema %>%
  ggplot(aes(x = frutos, y = proporcao)) +
  theme_bw() +
  geom_point(aes(pch = variedade2)) +
  geom_line(data = gema.pred2, aes(y = pred, lty = variedade2)) +
  xlab("Número de frutos") +
  ylab("Proporção de gemas florais")

# Gráficos do livro
par(mfrow=c(1,2))
hnp(m5, xlab="Quantis da distribuição meio-normal", ylab="Resíduos componentes do desvio", main="(a)", pch=16, cex=.7)
plot(proporcao ~ frutos, data=gema, pch=rep(16:18, each=5), cex=.9, ylim=c(0,.7),
     xlab="Número de frutos", ylab="Proporção de gemas florais", main="(b)")
legend(2,.72, c("Crispin","Cox","Golden Delicious"), inset=.01,
       bty="n", col=1, lwd=2, lty=1:3, pch=16:18, y.intersp=1.5, cex=.8)
lines(pred ~ frutos, data=subset(gema.pred, gema.pred$variedade=="Crispin"),
      col=1, lwd=2, lty=1)
lines(pred ~ frutos, data=subset(gema.pred, gema.pred$variedade=="Cox"),
      col=1, lwd=2, lty=2)
lines(pred ~ frutos, data=subset(gema.pred, gema.pred$variedade=="Golden Delicious"),
      col=1, lwd=2, lty=3)

par(mfrow=c(1,2))
hnp(m5.2, xlab="Quantis da distribuição meio-normal", ylab="Resíduos componentes do desvio", main="(a)", pch=16, cex=.7)
plot(proporcao ~ frutos, data=gema, pch=rep(16:18, each=5), cex=.9, ylim=c(0,.7),
     xlab="Número de frutos", ylab="Proporção de gemas florais", main="(b)")
legend(2,.72, c("Crispin","Cox","Golden Delicious"), inset=.01,
       bty="n", col=1, lwd=2, lty=c(1,2,2), pch=16:18, y.intersp=1.5, cex=.8)
lines(pred ~ frutos, data=subset(gema.pred2, gema.pred2$variedade=="1"),
      col=1, lwd=2, lty=1)
lines(pred ~ frutos, data=subset(gema.pred2, gema.pred2$variedade=="2"),
      col=1, lwd=2, lty=2)

########################################################################
## 7.5 Cultura de tecidos de macieiras
########################################################################
## Entrada dos dados
cult <- data.frame(tipo = gl(2, 90, labels = c("BAP","TDZ")),
                   nivel = gl(3, 30, labels = c("5.0","1.0","0.1")),
                   auxina = gl(3, 10, labels = c("NAA","IBA","2-4D")),
                   bloco = gl(10, 1),
                   regen = c(1,1,0,0,1,0,1,0,1,1,
                             0,1,1,1,1,1,0,1,1,1,
                             1,1,1,1,1,1,1,0,0,1,
                             0,0,0,0,0,0,0,0,0,0,
                             1,1,1,0,0,1,1,0,1,1,
                             1,0,1,1,0,1,1,1,1,1,
                             0,0,1,1,1,0,1,0,0,0,
                             0,0,0,1,1,1,1,0,1,0,
                             0,0,1,1,1,1,1,0,1,1,
                             1,1,1,1,1,0,1,1,1,1,
                             1,1,1,1,1,1,1,1,1,1,
                             1,0,1,1,1,1,1,1,1,1,
                             1,1,1,1,1,1,1,1,1,1,
                             1,1,1,1,1,1,1,1,1,1,
                             1,1,1,1,1,1,1,1,1,0,
                             1,1,1,1,1,1,1,0,1,1,
                             1,1,1,1,1,1,0,1,1,1,
                             0,0,1,0,1,1,1,1,1,1))

## Ajuste dos modelos

m1 <- glm(regen ~ 1, family = binomial, data = cult)
m2 <- glm(regen ~ bloco, family = binomial, data = cult)
m3 <- glm(regen ~ bloco + tipo, family = binomial, data = cult)
m4 <- glm(regen ~ bloco + nivel, family = binomial, data = cult)
m5 <- glm(regen ~ bloco + auxina, family = binomial, data = cult)
m6 <- glm(regen ~ bloco + tipo + nivel, family = binomial, data = cult)
m7 <- glm(regen ~ bloco + tipo + auxina, family = binomial, data = cult)
m8 <- glm(regen ~ bloco + nivel + auxina, family = binomial, data = cult)
m9 <- glm(regen ~ bloco + tipo * nivel, family = binomial, data = cult)
m10 <- glm(regen ~ bloco + tipo * auxina, family = binomial, data = cult)
m11 <- glm(regen ~ bloco + nivel * auxina, family = binomial, data = cult)
m12 <- glm(regen ~ bloco + tipo * nivel + auxina, family = binomial, data = cult)
m13 <- glm(regen ~ bloco + tipo * auxina + nivel, family = binomial, data = cult)
m14 <- glm(regen ~ bloco + nivel * auxina + tipo, family = binomial, data = cult)
m15 <- glm(regen ~ bloco + tipo * nivel + tipo * auxina, family = binomial, data = cult)
m16 <- glm(regen ~ bloco + tipo * nivel + nivel * auxina, family = binomial, data = cult)
m17 <- glm(regen ~ bloco + nivel * auxina + tipo * auxina, family = binomial, data = cult)
m18 <- glm(regen ~ bloco + (nivel + tipo + auxina)^2, family = binomial, data = cult)
m19 <- glm(regen ~ bloco + tipo * nivel * auxina, family = binomial, data = cult)

## g.l., desvios e X^2 (Tabela 7.16)
cv <- g.l. <- Desvios <- X2 <- NULL
for(i in 1:19) {
  modelo <- get(paste("m", i, sep = ""))
  cv[i] <- as.character(formula(modelo))[3]
  g.l.[i] <- modelo$df.residual
  Desvios[i] <- round(deviance(modelo), 2)
  X2[i] <- round(sum(resid(modelo, type = "pearson")^2),1)
}

tabela7.16 <- data.frame(cv, g.l., Desvios, X2)
tabela7.16

## Análises de desvios
anova(m19, test = "Chisq")
anova(m10, test = "Chisq")

## Análise de resíduos
require(hnp)
hnp(m10, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", pch = 16, cex = .7)

########################################################################
## 7.6 Toxicidade a dissulfeto de carbono gasoso
########################################################################
# Pacotes necessários
require(hnp)

# Entrada dos dados
# Collett (1991), página 109

y <- c(2,4,7,6,9,9,14,14,23,29,29,24,29,32,29,31)
m <- c(29,30,30,30,28,34,27,29,30,33,31,28,30,32,29,31)
dose <- rep(c(49.06,52.99,56.91,60.84,64.76,68.69,72.61,76.54), each = 2)

plot(y/m ~ dose, xlab = "Dose", ylab = "Proporção de insetos mortos")
resp <- cbind(y, m-y)

## Função de ligação logística
## regressão quadrática + desvios de regressão
modl <- glm(resp ~ poly(dose, 2) + factor(dose), family = binomial)
anova(modl, test="Chisq")

## Modelo nulo
mod1l <- glm(resp ~ 1, family = binomial)
pchisq(deviance(mod1l), df.residual(mod1l), lower.tail = FALSE)
print(sum(residuals(mod1l, 'pearson')^2)) # estatística X2

## Modelo linear
mod2l <- glm(resp ~ dose, family = binomial)
pchisq(deviance(mod2l), df.residual(mod2l), lower.tail = FALSE)
print(sum(residuals(mod2l, 'pearson')^2))

## Modelo quadrático
mod3l <- glm(resp ~ poly(dose, 2), family = binomial)
pchisq(deviance(mod3l), df.residual(mod3l), lower.tail = FALSE)
print(sum(residuals(mod3l, 'pearson')^2))

## Considerando dose como fator qualitativo
mod4l <- glm(resp ~ factor(dose), family = binomial)
pchisq(deviance(mod4l), df.residual(mod4l), lower.tail = FALSE)
print(sum(residuals(mod4l, 'pearson')^2))

## Testes para modelos encaixados
anova(mod1l, mod2l, mod3l, mod4l, test = "Chisq")

## Teste para função de ligação
LP2 <- predict(mod2l)^2
mod5l <- update(mod2l , . ~ . + LP2, family = binomial)
anova(mod5l, test="Chisq")

## Cálculo da concentração letal para o modelo logístico
coef_logistico <- coef(glm(resp ~ dose + I(dose^2), family = binomial))
round(polyroot(coef_logistico)[2], 4) ## CL50
coef_logistico2 <- coef_logistico
coef_logistico2[1] <- coef_logistico2[1] - qlogis(.9)
round(polyroot(coef_logistico2)[2], 4) ## CL90

## Função de ligação complemento log-log
modc <- glm(resp ~ poly(dose, 2) + factor(dose),
            family = binomial(link = "cloglog"))
anova(modc, test="Chisq")

## Modelo nulo
mod1c <- glm(resp ~ 1, family = binomial(link = "cloglog"))
pchisq(deviance(mod1c), df.residual(mod1c), lower.tail = FALSE)
print(sum(residuals(mod1c, 'pearson')^2))

## Modelo linear
mod2c <- glm(resp ~ dose, family = binomial(link = "cloglog"))
pchisq(deviance(mod2c), df.residual(mod2c), lower.tail = FALSE)
print(sum(residuals(mod2c, 'pearson')^2))

## Modelo quadŕatico
mod3c <- glm(resp ~ poly(dose,2), family=binomial(link = "cloglog"))
pchisq(deviance(mod3c), df.residual(mod3c), lower.tail = FALSE)
print(sum(residuals(mod3c, 'pearson')^2))

## Considerando dose como fator qualitativo
mod4c <- glm(resp ~ factor(dose), family=binomial(link = "cloglog"))
pchisq(deviance(mod4c), df.residual(mod4c), lower.tail = FALSE)
print(sum(residuals(mod4c, 'pearson')^2))

## Testando modelos encaixados
anova(mod1c, mod2c, mod3c, mod4c, test = "Chisq")

pdf("cs2_hnp_ambos.pdf", w = 12, h = 6)
par(mfrow = c(1,2))
hnp(mod3l, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", main = "(a)",
    pch = 16, cex = .7, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.6)
hnp(mod2c, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", main = "(b)",
    pch = 16, cex = .7, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.6)
dev.off()

## Cálculo da concentração letal para o modelo complemento log-log
dose.p(mod2c, p = c(.5, .9))

## Teste para função de ligação
LP2c <- predict(mod2c)^2
mod5c <- update(mod2c, . ~ . + LP2c, family = binomial(link = "cloglog"))
anova(mod5c, test="Chisq")

## Curvas ajustadas
pred <- expand.grid(dose = seq(49, 77, 0.1))
pred$plogit <- predict(mod3l, pred, type = "response")
pred$pcloglog <- predict(mod2c, pred, type = "response")

dados <- tibble(y = y,
                m = m,
                dose = dose)
names(pred)[2:3] <- c("logística","complemento log-log")

pdf("CS2_curvas_ajustadas.pdf", w = 6, h = 4)
dados %>%
  ggplot(aes(x = dose, y = y/m)) +
  theme_bw() +
  geom_point(alpha = .5) +
  geom_line(data = pred %>%
              pivot_longer(cols = 2:3,
                           names_to = "Função de ligação",
                           values_to = "pred"),
            aes(y = pred, lty = `Função de ligação`)) +
  xlab("Concentração de dissulfeto de carbono gasoso") +
  ylab("Proporção de insetos mortos")
dev.off()

########################################################################
## 7.7 Armazenamento de microorganismos
########################################################################
# Pacotes necessários
require(hnp)

# Entrada dos dados
bac <- data.frame(tempo = c(0,1,2,6,12),
                  contagem = c(31,26,19,15,20))

fit <- glm(contagem ~ log(tempo + .1), family = poisson, data = bac)

x <- seq(0, 12, length = 100)
pred <- predict(fit, data.frame(tempo = x), type = "response")

par(mfrow = c(1,2))
hnp(fit, xlab = "Quantis da distribuição meio-normal", ylab = "Resíduos componentes do desvio", main = "(a)", pch = 16, cex = .7)
plot(contagem ~ tempo, data = bac, xlab = "Tempo (meses)", ylab = "Contagens", 
     main = "(b)", pch = 16, cex = .8)
lines(x, pred)

## Ajustes
X2 <- function(obj) sum(resid(obj, type = "pearson")^2)

mod1 <- glm(contagem ~ 1, family = poisson, data = bac)
deviance(mod1)
pchisq(deviance(mod1), df.residual(mod1), lower.tail = FALSE)
X2(mod1)
pchisq(X2(mod1), df.residual(mod1), lower.tail = FALSE)

mod2 <- glm(contagem ~ log(tempo + .1), family = poisson, data = bac)
deviance(mod2)
pchisq(deviance(mod2), df.residual(mod2), lower.tail = FALSE)
X2(mod2)
pchisq(X2(mod2), df.residual(mod2), lower.tail = FALSE)

anova(mod2, test = "Chisq")

########################################################################
## 7.8 Número de brotos em um estudo de micropropagação de macieiras
########################################################################
micro.dat <- read.table("microprop.dat", header = TRUE)
for(j in 3:4) micro.dat[[j]] <- as.factor(micro.dat[[j]])

mod1 <- glm(Brotos ~ Meio * Hormonio, family = poisson, data = micro.dat)
1 - pchisq(deviance(mod1), df.residual(mod1))
anova(mod1, test="Chisq")

mod2 <- update(mod1, . ~ . + Erecip)
1 - pchisq(deviance(mod2), df.residual(mod2))
anova(mod1, mod2, test="Chisq")

require(hnp)
m1hnp <- hnp(mod1)
m2hnp <- hnp(mod2)

pdf("Fig-brotos-hnp.pdf", w = 6, h = 6)
plot(m2hnp, xlab = "Quantis da distribuição meio-normal", cex.lab = 1.5, cex.axis = 1.3,
     ylab = "Resíduos componentes do desvio", main = "", pch = 16, cex = .7)
dev.off()

########################################################################
## 7.9 Número de espécies de plantas
########################################################################
species <- read.table("especies.txt", header = TRUE)
levels(species$pH) <- c("Alto","Baixo","Médio")

spp <- with(species, split(Species, pH))
bio <- with(species, split(Biomass, pH))

species %>%
  ggplot(aes(x = Biomass, y = Species)) +
  theme_bw() +
  geom_point(aes(colour = pH)) +
  xlab("Biomassa") +
  ylab("Número de espécies")

model1 <- glm(Species ~ pH * Biomass, family = poisson, data = species)
deviance(model1)
pchisq(deviance(model1), df.residual(model1), lower.tail = FALSE)
model2 <- glm(Species ~ pH + Biomass, family = poisson, data = species)
deviance(model2)
pchisq(deviance(model2), df.residual(model2), lower.tail = FALSE)
model3 <- glm(Species ~ Biomass, family = poisson, data = species)
deviance(model3)
model4 <- glm(Species ~ 1, family = poisson, data = species)
deviance(model4)
anova(model4, model3, model2, model1, test = "Chisq")

anova(model1, test = "Chisq")
coef(update(model1, . ~ . - 1 - Biomass))

require(hnp)

par(mfrow = c(1,2), cex.axis = .8, cex.lab = .8)
hnp(model1, xlab = "Quantis da distribuição meio-normal", 
    ylab = "Resíduos componentes do desvio", main = "(a)", pch = 16, cex = .7)
with(species, plot(Biomass, Species, type = "n", xlab = "Biomassa",
                   ylab = "Número de espécies", main = "(b)"))
with(species, points(bio[[1]], spp[[1]], pch = 16, cex = .8))
with(species, points(bio[[2]], spp[[2]], pch = 17, cex = .8))
with(species, points(bio[[3]], spp[[3]], pch = 1, cex = .8))
legend("topright", c("Alto","Médio", "Baixo"), pch = c(16,1,17), bty = "n",
       col = 1, cex = .7, inset = .01, lty = c(1,2,4), lwd = 2, y.intersp = 1.5)
bio.x1 <- seq(0.087, 9.982, length = 30)
lines(bio.x1, predict(model1, data.frame(Biomass = bio.x1, pH = factor("Alto")), 
                      type = "response"), lwd = 2)
bio.x2 <- seq(0.176, 8.3, length = 30)
lines(bio.x2, predict(model1, data.frame(Biomass = bio.x2, pH = factor("Médio")), 
                      type = "response"), lty = 2, lwd = 2)
bio.x3 <- seq(0.05, 4.871, length = 30)
lines(bio.x3, predict(model1, data.frame(Biomass = bio.x3, pH = factor("Baixo")), 
                      type = "response"), lty = 4, lwd = 2)

########################################################################
## 7.10 Coleta de insetos em armadilhas adesivas
########################################################################
y <- c(246, 17, 458, 32)
armcor <- factor(c(1, 1, 2, 2))
sexo <- factor(c(1, 2, 1, 2))

count.dat <- data.frame(armcor, sexo, y)

# razão de chances observada
246*32/(17*458)

# ajuste do modelo loglinear
mod1 <- glm(y ~ armcor * sexo, family = poisson)
print(sum(residuals(mod1, 'pearson')^2))
anova(mod1, test = "Chisq")
summary(mod1)

# Note que este modelo reproduz os dados
# Também a razão de chances ajustada na escala log é 0.01098
exp(mod1$coef[4])

# A interação não é significativa, então não podemos rejeitar
# a hipótese de que a razão de chances é 1, isto é,
# cor de armadilha e sexo são independentes.

# Ajustando o modelo adequado -- mais simples
mod2 <- glm(y ~ armcor + sexo, family = poisson)
print(sum(residuals(mod1, 'pearson')^2))
anova(mod2, test = "Chisq")
1-pchisq(deviance(mod2), df.residual(mod2))
summary(mod2)

########################################################################
## 7.11 Pneumoconiose em mineiros de carvão
########################################################################
miners <- scan("Miners.dat", what = list(N = 0, M = 0, S = 0)) %>% as_tibble

miners <- miners %>%
  mutate(MS = M + S,
         Total = N + M + S,
         Ano_fator = factor(1:8),
         Ano = c(5.8,15,21.5,27.5,33.5,39.5,46,51.5))

miners_long <- miners %>%
  pivot_longer(1:3,
               names_to = "Categoria",
               values_to = "Contagens") %>%
  mutate(Categoria = factor(Categoria, levels = c("N","M","S")))

plot1 <- miners_long %>%
  ggplot(aes(x = Ano, y = Contagens/Total)) +
  theme_bw() +
  geom_point(aes(pch = Categoria)) +
  geom_line(aes(lty = Categoria)) +
  xlab("Número de anos de exposição") +
  ylab("Proporções de mineiros de carvão") +
  ggtitle("(a)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

plot2 <- miners %>%
  pivot_longer(2:3,
               names_to = "Categoria",
               values_to = "Contagens") %>%
  ggplot(aes(x = log(Ano), y = log(Contagens/N + .01))) +
  theme_bw() +
  geom_point(aes(pch = Categoria)) +
  geom_line(aes(lty = Categoria)) +
  xlab("log(Número de anos de exposição + 0.01)") +
  ylab("Logitos empíricos") +
  ggtitle("(b)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

library(gridExtra)
pdf("Fig_mineiros_exploratorio.pdf", w = 10, h = 4)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

# modelo nulo multinomial (reproduz totais marginais de ano e categoria)
mod1 <- glm(Contagens ~ Ano_fator + Categoria,
            family = poisson,
            data = miners_long) # mod indep.
summary(mod1)

mod2 <- glm(Contagens ~ Ano_fator + Categoria * log(Ano),
            family = poisson,
            data = miners_long)
anova(mod1, mod2, test = "Chisq")
summary(mod2)

library(hnp)
set.seed(1234)
hnp1 <- hnp(mod1)
hnp2 <- hnp(mod2)

pdf("mineiros_hnp.pdf", w = 12, h = 6)
par(mfrow = c(1,2))
plot(hnp1, xlab = "Quantis da distribuição meio-normal",
    ylab = "Resíduos componentes do desvio", main = "(a)", pch = 16, cex = .7)
plot(hnp2, xlab = "Quantis da distribuição meio-normal",
     ylab = "Resíduos componentes do desvio", main = "(b)", pch = 16, cex = .7)
dev.off()

## curvas preditas
miners_long$pred1 <- predict(mod1, type = "response")/miners_long$Total
miners_long$pred2 <- predict(mod2, type = "response")/miners_long$Total

plot1_pred <- miners_long %>%
  ggplot(aes(x = Ano, y = Contagens/Total)) +
  theme_bw() +
  geom_point(aes(pch = Categoria)) +
  geom_line(aes(y = pred1, lty = Categoria)) +
  xlab("Número de anos de exposição") +
  ylab("Proporções de mineiros de carvão") +
  ggtitle("(a)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

plot2_pred <- miners_long %>%
  ggplot(aes(x = Ano, y = Contagens/Total)) +
  theme_bw() +
  geom_point(aes(pch = Categoria)) +
  geom_line(aes(y = pred2, lty = Categoria)) +
  xlab("Número de anos de exposição") +
  ylab("Proporções de mineiros de carvão") +
  ggtitle("(b)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

library(gridExtra)
pdf("Fig_mineiros_predito.pdf", w = 10, h = 4)
grid.arrange(plot1_pred, plot2_pred, ncol = 2)
dev.off()