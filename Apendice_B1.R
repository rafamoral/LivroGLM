# Entrando com os dados
d <- c(0.0, 2.6, 3.8, 5.1, 7.7, 10.2) # doses
m <- c(49, 50, 48, 46, 49, 50) # números de insetos que receberam as doses
y <- c(0, 6, 16, 24, 42, 44) # números de insetos mortos

# Modelo binomial -- Y ~ Bin(m, pi) -- com preditor linear
# eta = logit(pi) = beta1 + beta2*dose

# Primeira iteração - beta^{(1)}
y[1] <- .1              # Para prevenir problemas numéricos
mu <- y                 # mu = y
eta <- log(y/(m-y))     # eta = logit(mu) = logit(y)
vy <- y*(m-y)/m         # V(mu) = (mu/m)*(m - mu) = (y/m)*(m - y)
glinha.y <- m/(y*(m-y)) # g'(mu) = g'(y)
w <- 1/(vy*glinha.y^2)  # w = 1/(V(mu) * g'(mu)^2) = 1/(V(y) * g'(y)^2)
W <- diag(w)            # W = diag(w)
X <- model.matrix(~ d)  # matriz de delineamento

beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% eta # beta^{(2)}

# Da segunda iteração em diante:
i <- 1
crit <- 1
solucoes <- data.frame(beta)

while(crit > 1e-16 & i < 20) {
  eta <- beta[1] + d*beta[2]                                # preditor linear
  mu <- m*exp(eta)/(1 + exp(eta))                           # mu
  glinha <- m/(mu*(m-mu))                                   # g'(mu)
  z <- eta + (y - mu)*glinha                                # z
  vmu <- mu*(m-mu)/m                                        # V(mu)
  w <- 1/(vmu * glinha^2)                                   # w
  W <- diag(w)                                              # W = diag(w)
  beta.novo <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z # beta da próxima iteração
  i <- i + 1
  crit <- sum(((beta.novo - beta)/beta)^2)                  # critério de convergência
  cat("iteracao", i, ": conv =", crit, "\n")
  solucoes <- data.frame(solucoes, beta.novo)
  beta <- beta.novo                                         # atualizando beta
}

names(solucoes) <- paste("iter", 1:i, sep = "")

solucoes # valores dos betas para cada iteração

# Utilizando o comando glm()
y[1] <- 0
glm(cbind(y, m-y) ~ d, family = binomial) # a função de ligação logit é a
                                          # padrão quando se usa family = binomial
