normalize <- function(x){
  x/sqrt(sum(x^2))
}

V <- rPKBD_Saw(1000, 0.99, c(0,0,1))
plot(density(V[,3]))
V <- rPKBD_ACG(1000, 0.99, c(0,0,1))
lines(density(V[,3]), col = 2)

V <- rPKBD_Saw(1000, 0.4, c(0,0,1))
plot(density(V[,3]))
V <- rPKBD_ACG(1000, 0.4, c(0,0,1))
lines(density(V[,3]), col = 2)

V <- rPKBD_Saw(1000, 0.1, c(0,0,0,1))
plot(density(V[,4]))
V <- rPKBD_ACG(1000, 0.1, c(0,0,0,1))
lines(density(V[,4]), col = 2)

V <- rPKBD_Saw(1000, 0.7, normalize(1:5))
plot(density(V[,3]))
V <- rPKBD_ACG(1000, 0.7, normalize(1:5))
lines(density(V[,3]), col = 2)

V <- rPKBD_Saw(1000, 0.3, normalize(rep(1,15)))
plot(density(V[,6]))
V <- rPKBD_ACG(1000, 0.3, normalize(rep(1,15)))
lines(density(V[,6]), col = 2)

V <- rPKBD_Saw(1000, 0.99, normalize(1:200))
plot(density(V[,6]))
V <- rPKBD_ACG(1000, 0.99, normalize(1:200))
lines(density(V[,6]), col = 2)

V <- rPKBD_Saw(1000, 0.99, normalize(1:2))
plot(density(V[,2]))
V <- rPKBD_ACG(1000, 0.99, normalize(1:2))
lines(density(V[,2]), col = 2)

V <- rPKBD_Saw(1000, 0.99, c(0,0,1))
library(rgl)
rgl.open()
view3d(zoom =0.5)
rgl.bg(color = "white")
rgl.points(V[,1],V[,2],V[,3], ylim=c(-1,1),col = "red" ,xlim=c(-1,1), zlim = c(-1,1), xlab = "x", ylab = "y", zlab = "z", alpha = 0.8)
spheres3d(x = 0, y = 0, z = 0, radius = 0.98, col = "green", alpha = 0.6, back = "lines")
rgl.lines(c(-1.5,1.5), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(-1.5,1.5), c(0, 0), color = "blue")
rgl.lines(c(0, 0), c(0, 0), c(-1.5,1.5), color = "red")

V <- rPKBD_ACG(1000, 0.99, c(0,0,1))
rgl.open()
view3d(zoom =0.5)
rgl.bg(color = "white")
rgl.points(V[,1],V[,2],V[,3], ylim=c(-1,1),col = "red" ,xlim=c(-1,1), zlim = c(-1,1), xlab = "x", ylab = "y", zlab = "z", alpha = 0.8)
spheres3d(x = 0, y = 0, z = 0, radius = 0.98, col = "green", alpha = 0.6, back = "lines")
rgl.lines(c(-1.5,1.5), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(-1.5,1.5), c(0, 0), color = "blue")
rgl.lines(c(0, 0), c(0, 0), c(-1.5,1.5), color = "red")