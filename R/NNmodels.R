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