
library(survival)
library(frailtyHL)
library(frailtyEM)

gf = GaFrailtyMM(Surv(time, status) ~ age + sex + cluster(id), data=kidney, Maxiter = 1000, convergence = 1e-06)
summary.GaFrailtyMM(gf)

# gf$print_k

# frailtyHL
kidney_g02<-frailtyHL(Surv(time,status)~sex+age+(1|id),kidney,RandDist="Gamma",mord=0,dord=2)

# coxph
coxph(Surv(time, status) ~ age + sex + frailty.gamma(id, eps=1e-6, method="em", sparse=0),
             outer.max=1000, iter.max=10000,
             data=kidney)

# frailtyEM
emfrail(formula = Surv(time, status) ~ sex + age + cluster(id),
        data = kidney, distribution = emfrail_dist(dist = "gamma"))

