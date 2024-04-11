# LTNMM

set.seed(1)
y <- LTN_sample(5, 7, 4, 100)
result <- LTNMM(y ~ 1, a = 5)
summary(result)


# ZIGPMM

set.seed(1)
m <- 2
phi0 <- 0.1
la <- rep(9, m)
th <- rep(0.7, m)

y <- ZIGP_sample(phi0, la, th, 600)

la <- rep(1, m)
th <- rep(0.1, m)

result <- ZIGPMM(y, phi0, la, th)
summary(result)

# CZIGPMM

set.seed(1)
n <- 600
m <- 2
phi0 <- 0.1
phi <- c(0.2, 0.5)
la <- c(2, 7)
th <- c(0.2, 0.3)
y <- CZIGP_sample(phi0, phi, la, th, n)

phi0 <- 0.1
phi <- rep(0.1, m)
la <- rep(1, m)
th <- rep(0.1, m)

result <- CZIGPMM(y, phi0, phi, la, th)
summary(result)

# COXMM
#

#
# # IC2MM
#
# sample <- function(la, n) {
#   u <- runif(n, 0, 1)
#   w <- runif(n, 0, 1)
#   U <- -log(u)/la
#   W <- -log(w)/la
#   L <- pmin(U, W)
#   R <- pmax(U, W)
#   return(list(L = L, R = R))
# }
# n <- 10
# y <- sample(0.5, n)
# L <- y$L
# R <- y$R
# data <- data.frame(L, R)
#
# result <- IC2MM(Surv(L, R, type = "interval2") ~ 1, data, control = IC2Control(Pdigits = 3))
# summary(result)
# plot(result)
#
# # 2
# result <- IC2MM(Surv(left, right, type = "interval2") ~ treatment, bcos,
#                 IC2Control(Pdigits = 3))
# summary(result)
# plot(result, col = c("red", "blue"))
#
# # 2
#
# phi0 <- 0.5
# la <- rep(5, 2)
# th <- rep(0.5, 2)
#
# result <- ZIGPMM(vijc, phi0, la, th)
# summary(result)
#
# # 3
# result <- ZIGPMM(cadi, phi0, la, th)
# summary(result)

library(survival)

# CoxMM

# 1
set.seed(1)
n <- 200
q <- 3
be <- matrix(c(-0.5, 1, 2), ncol = 1)
t <- c(0.1, 0.5)

da <- COXsample(n, be = be, la = 2)
data <- da$x
time <- da$t
status <- da$I

data <- data.frame(data[, 1], data[, 2], data[, 3], time = time, status)

result <- CoxMM(Surv(time, status) ~ X1 + X2 + X3, data, beat = be)
summary(result)
plot(result)

# 2
result <- CoxMM(Surv(time, status) ~ age + sex, lung)
summary(result)
plot(result)

# GaFrailtyMM

result <- GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data = kidney)
summary(result)
plot(result)

