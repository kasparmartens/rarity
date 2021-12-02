import torch
import torch.nn as nn
import numpy as np
from torch.distributions.normal import Normal
from torch.distributions.relaxed_bernoulli import LogitRelaxedBernoulli

from rarity.helpers import summarise_binary_profiles


class RelaxedBernoulliAutoEncoder(nn.Module):
    """
    Autoencoder with a binary latent space, where dim(z) == dim(data),
    fitted via continuous relaxations of binary random variables
    """
    def __init__(self, data_dim, optimise_params=False):
        super().__init__()

        self.optimise_params = optimise_params

        self._init_decoder_params()

        self.encoder = nn.Sequential(
            nn.Linear(data_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, data_dim)
        )

    def _init_decoder_params(self):
        if self.optimise_params:
            self.mu1_logit = nn.Parameter(torch.zeros(1))
            self.mu2_logit = nn.Parameter(torch.zeros(1))
            self.sigma1_logit = nn.Parameter(torch.zeros(1))
            self.sigma2_logit = nn.Parameter(torch.zeros(1))
        else:
            # option to keep params fixed
            self.mu1 = 0.05
            self.mu2 = 0.5
            self.sigma1 = 0.03
            self.sigma2 = 0.2

    def decode(self, z):
        if self.optimise_params:
            # mu1 restricted to [0, 0.1]
            mu1 = 0.1 * torch.sigmoid(self.mu1_logit)
            # mu2 restricted to [0.4, 0.6]
            mu2 = 0.4 + 0.2 * torch.sigmoid(self.mu2_logit)
            # sigma1 restricted to [0.01, 0.06]
            sigma1 = 0.01 + 0.05 * torch.sigmoid(self.sigma1_logit)
            # sigma2 restricted to [0.1, 0.3]
            sigma2 = 0.1 + 0.2 * torch.sigmoid(self.sigma2_logit)
        else:
            mu1 = self.mu1
            mu2 = self.mu2
            sigma1 = self.sigma1
            sigma2 = self.sigma2
        mu = (1.0 - z) * mu1 + z * mu2
        sigma = (1.0 - z) * sigma1 + z * sigma2
        return mu, sigma

    def forward(self, Y, temperature):
        # encode
        logits = self.encoder(Y)
        # sample z
        p = LogitRelaxedBernoulli(temperature, logits=logits)
        z_sample = torch.sigmoid(p.rsample())
        # decoding
        mu, sigma = self.decode(z_sample)
        return mu, sigma, logits

    def loglik(self, Y, mu, sigma):
        return Normal(loc=mu, scale=sigma).log_prob(Y)


def fit_Rarity(Y, n_iter=20000, batch_size=1024, optimise_params=False, verbose=True):
    """
    Procedure for fitting the Rarity clustering model.
    :param Y: The observed expression matrix with cells in rows and genes in columns. We assume the intensities have been scaled to [0, 1].
    :param n_iter: The number of iterations.
    :param batch_size: The number of cells in a batch.
    :param optimise_params: Indicates whether the likelihood parameters are fixed or optimised
    :param verbose:
    :return:
    """
    Y = torch.Tensor(Y)
    N = Y.shape[0]

    temperature_grid = torch.linspace(4.0, 0.2, steps=n_iter)
    m = RelaxedBernoulliAutoEncoder(data_dim=Y.shape[1], optimise_params=optimise_params)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3)

    for i in range(n_iter):
        subsample_idx = np.random.permutation(N)[:batch_size]
        Y_sub = Y[subsample_idx, :]
        mu, sigma, logits = m.forward(Y_sub, temperature_grid[i])
        loglik = m.loglik(Y_sub, mu, sigma)
        loss = - loglik.sum()
        opt.zero_grad()
        loss.backward()
        opt.step()
        if verbose and (i % 1000 == 0):
            print(f"iter {i:5d} loss {loss.item():1.2f}")

    with torch.no_grad():
        _, _, logits = m.forward(Y, temperature_grid[-1])
        z_mean = torch.sigmoid(logits).numpy()
        z_mean_binary = 1 * (z_mean > 0.5)

    counts, unique_profiles, cluster_allocations = summarise_binary_profiles(z_mean_binary)

    return {
        "z_mean": z_mean,
        "z_binary": z_mean_binary,
        "counts": counts,
        "unique_profiles": unique_profiles,
        "cluster_allocations": 1 + cluster_allocations
    }
