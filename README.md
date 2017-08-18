# Bayesian-Probing
This repository consists of documents necessary to implement the probing methods to determine sparsity of a matrix S.
The methods include Griewank and Mitev's single and batch probing methods, and the batch probing method we suggested.

The single probing method can be implemented and the number of probes required to determine the sparsity can be found using find_probe.m
find_probe.m calls functions single_probe.m to determine the support of the single probe and update_probability.m to update the probability matrix P iteratively.

Batch_probing_methods.m implements the batch probing mehods. The method suggested by Griewank and Mitev calls functions single_probe.m and update_probability_bundle.m iteratively till the complete sparsity pattern is determined. Our method calls functions bundle_probe2.m and update_probaability_bundle.m iteratively till the same convergence of P and S.
