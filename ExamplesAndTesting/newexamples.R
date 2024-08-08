library(flexmix)
library(circlus)

mix <- rbind(rPKBD_ACG(300, 0.95, c(1,0,0)), rPKBD_ACG(300, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = PKBDNN_clust())
