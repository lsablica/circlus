import torch
import torch.nn as nn
import torch.optim as optim
from scipy.special import ive
import numpy as np

class VMFModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(VMFModel, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.fc_mu_kappa = nn.Linear(hidden_dim, output_dim)
        
    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        mu_kappa = self.fc_mu_kappa(x)
        mu = mu_kappa / mu_kappa.norm(dim=-1, keepdim=True)  # Ensure mu is a unit vector
        kappa = mu_kappa.norm(dim=-1, keepdim=True) / (mu_kappa.norm(dim=-1, keepdim=True) + 1)  # Map norm to (0, 1)
        return mu, kappa

def vmf_nll(mu, kappa, x, responsibilities):
    d = mu.shape[-1]
    kappa_np = kappa.detach().cpu().numpy()  # Convert to numpy for Bessel function
    C_d = (kappa_np ** (d/2 - 1)) / ((2 * np.pi) ** (d/2) * ive(d/2 - 1, kappa_np))  # Normalizing constant
    C_d = torch.tensor(C_d, dtype=torch.float32, device=kappa.device)
    log_prob = torch.log(C_d) + kappa * (mu * x).sum(dim=-1)
    weighted_log_prob = responsibilities * log_prob
    return -weighted_log_prob.mean()

def e_step(X, pi, mu, kappa):
    responsibilities = []
    for i in range(len(pi)):
        C_d = (kappa[i] ** (X.shape[-1]/2 - 1)) / ((2 * np.pi) ** (X.shape[-1]/2) * ive(X.shape[-1]/2 - 1, kappa[i]))
        log_prob = torch.log(C_d) + kappa[i] * (mu[i] * X).sum(dim=-1)
        responsibilities.append(pi[i] * torch.exp(log_prob))
    responsibilities = torch.stack(responsibilities, dim=1)
    responsibilities = responsibilities / responsibilities.sum(dim=1, keepdim=True)
    return responsibilities

def m_step(X, responsibilities, model, optimizer):
    optimizer.zero_grad()
    mu, kappa = model(X)
    loss = vmf_nll(mu, kappa, X, responsibilities)
    loss.backward()
    optimizer.step()
    return mu, kappa

# Example usage
input_dim = 10  # Number of covariates
hidden_dim = 50
output_dim = 3  # Dimension of the sphere, e.g., 3 for 3D sphere
num_components = 3

model = VMFModel(input_dim, hidden_dim, output_dim)
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Initialize parameters
pi = torch.ones(num_components) / num_components
mu = torch.randn(num_components, output_dim)
mu = mu / mu.norm(dim=1, keepdim=True)
kappa = torch.ones(num_components)

# Training loop for EM algorithm
num_epochs = 100
for epoch in range(num_epochs):
    model.train()
    # Assume X is the input covariates and Y is the observed directions
    X = torch.randn(100, input_dim)
    Y = torch.randn(100, output_dim)
    Y = Y / Y.norm(dim=-1, keepdim=True)  # Ensure Y is on the unit sphere
    
    # E-step
    responsibilities = e_step(Y, pi, mu, kappa)
    
    # M-step
    mu, kappa = m_step(Y, responsibilities, model, optimizer)
    
    if epoch % 10 == 0:
        print(f'Epoch {epoch}, mu: {mu}, kappa: {kappa}')

# After training, you can use the learned parameters (mu, kappa) to estimate the von Mises-Fisher distribution.