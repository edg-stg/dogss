Hi! This is the github repository of the R package "dogss", which is short for DOuble-Group Sparse Spike-and-slab, a Bayesian feature selection method. It takes advantage of known group structure of the features, where sparsity is assumed on the between-group and within-group levels. Thus it is very similar to the idea of the "sparse-group lasso". Implementation is much faster than your typical Bayesian method, since we employed a deterministic algorithm called expectation propagation instead of the stochastic and slow Gibbs sampling.

For more details, see: https://arxiv.org/abs/1809.09367

We strongly recommend to check out the provided notebook(s) which give an idea what the package is for. Furthermore they give instructions how to reproduce the figures from the paper linked above: https://github.com/edgarst/dogss/blob/master/notebooks/signal_recovery/signal_recovery.md
