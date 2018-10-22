# dogss

Hi! This is the github repository of the R package "dogss", which is short for **do**uble-**g**roup **s**parse **s**pike-and-slab, a Bayesian feature selection method. It takes advantage of known group structure of the features, where sparsity is assumed on the between-group and within-group levels. Thus it is very similar to the idea of the [sparse-group lasso](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2012.681250#.W82hp59BrmE), but in contrast it is based on Bayesian parameter inference. Our implementation is much faster than your typical Bayesian method, since we employed a deterministic algorithm called expectation propagation instead of the stochastic and slow Gibbs sampling.

## Publication

For more details, see the corresponding publication [Sparse-Group Bayesian Feature Selection Using Expectation Propagation for Signal Recovery and Network Reconstruction](https://arxiv.org/abs/1809.09367).

## Notebooks

We strongly recommend to check out the provided notebook(s) which give an idea what the package is for. Furthermore they give instructions how to reproduce the figures from the paper linked above:

[Signal recovery with dogss (Figures 2-8)](notebooks/signal_recovery/signal_recovery.md)
