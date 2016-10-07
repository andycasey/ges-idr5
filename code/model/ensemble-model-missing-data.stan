
// Single-parameter ensemble model with correlation coefficients between nodes

data {
    int<lower=1> N; // number of nodes (estimators)
    int<lower=1> C; // number of calibrators (benchmarks)
    int<lower=1> V; // maximum number of visits to any calibrator
    int<lower=1> visits[C]; // number of visits to each calibrator
    
    int<lower=0> N_missing; // total number of missing data points.
    int is_missing[C, N, V]; // whether a node estimate is missing (>0 = missing)

    matrix[N, V] estimates[C]; // estimates of the stellar property
    matrix[V, C] spectrum_isnr; // inverse-snr of the spectrum (1/SNR)  

    real mu_calibrator[C]; // the non-spectroscopic or high fidelity calibrator value
    real sigma_calibrator[C]; // the 1\sigma uncertainty on the calibrator value

    real lower_alpha_sq; // lower bound on alpha parameter for f(snr)
    real upper_alpha_sq; // upper bound on alpha parameter for f(snr)

    real lower_bound; // lower univate bound on missing data points
    real upper_bound; // upper univate bound on missing data points
}
transformed data {
    int TCV; // total number of calibrator visits
    TCV = sum(visits);
}

parameters {
    real truths[C]; // God's word: true values of the calibrators.
    real biases[N]; // biases in the individual nodes
    real<lower=lower_bound, upper=upper_bound> missing_estimates[N_missing];

    // Cholesky factor of a correlation matrix
    cholesky_factor_corr[N] L_corr;
    
    // \alpha:  parameter to model uncertainty as a function of SNR s.t.
    //          \_variance_{random} = \alpha^2 * spectrum_isnr
    vector<lower=lower_alpha_sq, upper=upper_alpha_sq>[N] alpha_sq; 
    vector<lower=0, upper=upper_alpha_sq>[N] variance_sys;
}

transformed parameters {
    cov_matrix[N] Sigma[TCV]; // covariance matrix
    matrix[N, V] full_rank_estimates[C]; // array containing known (data) and 
                                         // unknown (parameter) estimates 
    {
        int mi; // an index for unknown data points
        int vi; // an index for total visits
        mi = 1;
        vi = 1;

        for (c in 1:C) { // For each calibrator
            for (v in 1:visits[c]) { // For each visit of this calibrator
                vector[N] sigma;

                // Build the diagonal of the covariance matrix
                for (n in 1:N) {
                    sigma[n] = sqrt(alpha_sq[n] * spectrum_isnr[v, c] +
                                    variance_sys[n]);
                
                    if (is_missing[c, n, v]) {
                        full_rank_estimates[c, n, v] = missing_estimates[mi] - biases[n];
                        mi = mi + 1;
                    }
                    else {
                        full_rank_estimates[c, n, v] = estimates[c, n, v] - biases[n];
                    }
                }
                // Construct the covariance matrix from Cholesky factorisation
                Sigma[vi] = diag_pre_multiply(sigma, L_corr) 
                          * diag_pre_multiply(sigma, L_corr)';
                vi = vi + 1;
            }
        }
    }
}

model {
    mu_calibrator ~ normal(to_vector(truths), sigma_calibrator); 

    L_corr ~ lkj_corr_cholesky(2);

    {
        int vi;
        vi = 1;

        for (c in 1:C) { // For each calibrator
            for (v in 1:visits[c]) { // For each visit of this calibrator
                full_rank_estimates[c, :, v] ~ multi_normal(rep_vector(truths[c], N)',
                                                            Sigma[vi]);

                vi = vi + 1;
            }
        }
    }

}

generated quantities {
    // For convenience.
    //real<lower=0> intrinsic_sigma;
    //vector<lower=0>[N_estimators] estimator_sys_sigma;
    //vector<lower=0>[N_estimators] estimator_rand_sigma;

    //intrinsic_sigma = pow(intrinsic_var, 0.5);
}
