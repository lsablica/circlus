library(torch)

convnet <- nn_module(
  "convnet",
  
  initialize = function() {
    
    # nn_conv2d(in_channels, out_channels, kernel_size)
    self$conv1 <- nn_conv2d(1, 16, 3)
    self$conv2 <- nn_conv2d(16, 32, 3)
    self$conv3 <- nn_conv2d(32, 64, 3)
    
    self$output <- nn_linear(2304, 3)
    
  },
  
  forward = function(x) {
    
    x %>% 
      self$conv1() %>% 
      nnf_relu() %>%
      nnf_max_pool2d(2) %>%
      self$conv2() %>% 
      nnf_relu() %>%
      nnf_max_pool2d(2) %>%
      self$conv3() %>% 
      nnf_relu() %>%
      nnf_max_pool2d(2) %>%
      torch_flatten(start_dim = 2) %>%
      self$output()
    
  }
)

model <- convnet()

img <- torch_randn(1, 1, 64, 64)

model(img)



train_dl <- dataloader(train_ds,
                       batch_size = 128,
                       shuffle = TRUE
)
valid_dl <- dataloader(valid_ds, batch_size = 128)



fitted <- convnet %>%
  setup(
    loss = nn_cross_entropy_loss(),
    optimizer = optim_adam,
    metrics = list(
      luz_metric_accuracy()
    )
  ) %>%
  fit(train_dl,
      epochs = 50,
      valid_data = valid_dl,
      verbose = TRUE
)



PKBD <- nn_module(
  "PKBD",
  initialize = function(input_dim, output_dim) {
    self$fc = nn_linear(input_dim, output_dim, bias = FALSE)
    self$output_dim = output_dim
  },
  forward = function(x) {
    mu = self$fc(x)
    normm = mu$norm(dim=-1, keepdim=TRUE)
    rho = normm / (1 + normm)
    mu = mu / normm
    list(mu = mu, rho = rho)
  }
)

pkbd_log_likelihood <- function(mu, rho, Y, W){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood = (d/2)*term2$log() - term1
  result = (neg_log_likelihood * W)$sum()
  return(result)
}


# Example usage
input_dim = 1
output_dim = 3
EPOCHS = 100
LR = 0.5
batch_size = 32



library(Rcpp) 
sourceCpp("~/Documents/GitHub/PKBD---code/src/rpkbd.cpp")
Y = rPKBD_ACG(1000, 0.95, c(1,0,0))
Y = torch_tensor(Y, dtype = torch_float32())
X = torch_ones(Y$shape[1], input_dim, dtype = torch_float32())
W = torch_ones(Y$shape[1], 1, dtype = torch_float32())/1000




model = PKBD(input_dim, output_dim)
optimizer = optim_adam(model$parameters, lr = LR)

for(epoch in seq_len(EPOCHS)){
  optimizer$zero_grad()
  res = model(X)
  print(res$rho[1]$item())
  loss = pkbd_log_likelihood(res$mu, res$rho, Y, W)
  print(loss$item())
  loss$backward()
  optimizer$step()
}


model = PKBD(input_dim, output_dim)
optimizer = optim_lbfgs(model$parameters, lr = LR, max_iter = 20, line_search_fn = "strong_wolfe")


calc_loss <- function() {
  optimizer$zero_grad()
  res = model(X)
  print(res$rho[1]$item())
  loss = pkbd_log_likelihood(res$mu, res$rho, Y, W)
  print(loss$item())
  loss$backward()
  loss
}
for(epoch in seq_len(2)){
  cat("\nIteration: ", epoch, "\n")
  optimizer$step(calc_loss)
}



# Example usage
input_dim = 2
output_dim = 3
EPOCHS = 2
LR = 0.5
batch_size = 32

Y = rbind(rPKBD_ACG(1000, 0.95, c(1,0,0)), rPKBD_ACG(2000, 0.7, c(0,0,1)))
Y = torch_tensor(Y, dtype = torch_float32())
X = torch_tensor(cbind(rep(1,3000), c(rep(1,1000),rep(0,2000))),  dtype = torch_float32()) 
W = torch_ones(Y$shape[1], 1, dtype = torch_float32())/3000




model = PKBD(input_dim, output_dim)
optimizer = optim_lbfgs(model$parameters, lr = LR, max_iter = 200, line_search_fn = "strong_wolfe")


calc_loss <- function() {
  optimizer$zero_grad()
  res = model(X)
  print(res$rho[1]$item())
  loss = pkbd_log_likelihood(res$mu, res$rho, Y, W)
  print(loss$item())
  loss$backward()
  loss
}
for(epoch in seq_len(1)){
  cat("\nIteration: ", epoch, "\n")
  optimizer$step(calc_loss)
}
res = model(X)
res$mu[200]; res$mu[2000]; res$rho[200]; res$rho[2000]




















res = model(X)
mu = res$mu
rho = res$rho
pkbd_log_likelihood(mu, rho, Y, W)




def pkbd_log_likelihood(mu, rho, x, W):
  d = x.shape[-1]
term1 = (1 - rho ** 2).log()  # log(1 - rho^2)
term2 = 1 + rho ** 2 - 2 * rho * ((mu.unsqueeze(1) @ Y.unsqueeze(2)).squeeze(-1))  # (1 + rho^2 - 2*rho*mu^T x)^(d/2)
log_likelihood = (d/2)*term2.log() - term1
#print(log_likelihood.shape, W.shape)
return (log_likelihood * W).sum()

# Example usage
input_dim = 1
output_dim = 3
EPOCHS = 100
LR = 0.5
batch_size = 32


Y = np.load('Y.npy')
Y = torch.tensor(Y, dtype=torch.float32)
Y = Y[:1000, :]
Y.shape

X = torch.ones(Y.shape[0], input_dim, dtype=torch.float32)
X = X[:1000, :]
X

W = torch.ones(Y.shape[0], input_dim, dtype=torch.float32)/1000
W.shape
W


model = PKBD(input_dim, output_dim)
optimizer = optim.Adam(model.parameters(), lr=LR)

%%time
model.train()
for epoch in range(EPOCHS):
  optimizer.zero_grad()
mu, rho = model(X)
print(rho[1].item())
loss = pkbd_log_likelihood(mu, rho, Y, W)
print(loss.item())
loss.backward()
optimizer.step()


model.fc.weight
