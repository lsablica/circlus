library(flexmix)
library(circlus)

mix <- rbind(rPKBD_ACG(300, 0.95, c(1,0,0)), rPKBD_ACG(300, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = circlus::PKBDNN_clust())
m1

df = readRDS("~/Documents/GitHub/PKBD---code/ExamplesAndTesting/df_final.RDS")
df = data.frame(df)
head(df)
names(df)[c(279,280)]


GTE <-  t(sapply(df[,279], function(x) x))
OAI <-  t(sapply(df[,280], function(x) x))
pages = matrix(df$pages, ncol = 1)

PKBD_abstract_6 <- flexmix(OAI ~ 1, k = 6, model = circlus::PKBD_clust())
PKBD_abstract_6
PKBD_abstract_6@logLik
sapply(PKBD_abstract_6@components, function(x) x[[1]]@parameters$rho)

PKBDNN_abstract_6 <- flexmix(OAI ~ 1, k = 6, model = circlus::PKBDNN_clust(LR= 0.01),
                             control = list(verbose = 1))
PKBDNN_abstract_6
PKBDNN_abstract_6@logLik

PKBDNN_abstract_8 <- flexmix(OAI ~ 1, k = 8, model = circlus::PKBDNN_clust(LR = 0.05),
                              control = list(verbose = 1))

PKBDNN_abstract_6b <- flexmix(OAI ~ 1 + pages, k = 6, model = circlus::PKBDNN_clust(LR = 0.005, max_iter = 20),
                              control = list(verbose = 1, nrep = 10))


