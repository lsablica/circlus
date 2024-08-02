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
    self$fc = nn_linear(input_dim, output_dim)
    self$output_dim = output_dim
  },
  forward = function(x) {
    mu = self.fc(x)
    normm = mu$norm(dim=-1, keepdim=True)
    rho = normm / (1 + normm)
    mu = mu / normm
    list(mu, rho)
  }
)


# Example usage
input_dim = 1
output_dim = 3
EPOCHS = 100
LR = 0.5
batch_size = 32

model = PKBD(input_dim, output_dim)


library(Rcpp) 
sourceCpp("~/Documents/GitHub/PKBD---code/src/rpkbd.cpp")
Y = rPKBD_ACG(1000, 0.95, c(1,0,0))
Y = torch_tensor(Y, dtype = torch_float32())




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

