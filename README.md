The code in this repository can be used to replicate results presented in ``Review and Demonstration of a Mixture Representation for Simulation from Densities Involving Sums of Powers.'' 

The files `powerregfixlamq_naive.stan`, `powerregfixlamq_centered.stan`, and `powerregfixlamq_noncentered.stan` provide `STAN` code for simulating a random vector $\boldsymbol z_2$ from log posteriors of the form 

$$
\frac{1}{2\sigma^2}\left\|\boldsymbol y - \boldsymbol X \boldsymbol z_2\right\|^2_2 + \lambda \left\|\boldsymbol z_2\right\|^q_q,
$$
for fixed $\sigma$, $\boldsymbol y$, $\boldsymbol X$, $\lambda$, and $q$ which refer to the noise standard deviation $\sigma > 0$, an $m\times 1$ response $\boldsymbol y$, an $m\times n_2$ matrix of covariates $\boldsymbol X$, a scale parameter $\lambda > 0$, and a shape parameter $0 < q < 2$.

The differences between the files are described in the corresponding note.

The file `simple_examination.R` simulates from log posteriors of the above form for two different datasets, which are provided in the `Data` folder in this directory. 

The file `summarize_simple_examination.R` reproduces plots provided in the note based on the results produced by the file `simple_examination.R`.
