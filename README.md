# dogss

Hi! This is the github repository of the R package "dogss", which is short for **do**uble-**g**roup **s**parse **s**pike-and-slab, a Bayesian feature selection method. It takes advantage of known group structure of the features, where sparsity is assumed on the between-group and within-group levels. Thus it is very similar to the idea of the sparse-group lasso, but in contrast it is based on Bayesian parameter inference. Our implementation is much faster than your typical Bayesian method, since we employed a deterministic algorithm called expectation propagation instead of the stochastic and slow Gibbs sampling.

## Publication

For more details, see the corresponding publication [Sparse-Group Bayesian Feature Selection Using Expectation Propagation for Signal Recovery and Network Reconstruction](https://arxiv.org/abs/1809.09367).

## Notebooks

We strongly recommend to check out the provided notebooks which give an idea what the package is for. Furthermore they give instructions how to reproduce the figures from the paper linked above:

* [Signal recovery with dogss (Figures 2-8)](notebooks/signal_recovery/signal_recovery.md)

* [Network reconstruction with dogss: Simulations (Figures 9-13)](notebooks/networks/networks.md)

* [Network reconstruction with dogss: DREAM5 (Figures 14 and 15)](notebooks/networks/networks_dream5.md)
