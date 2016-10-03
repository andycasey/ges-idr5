
// Mean ensemble model with correlation coefficients between nodes

functions {
    int num_pairwise(int N) {
        // Return the number of pair-wise comparisons for N items:
        // N_pairwise = N!/(2(N - 2)!)

        int N_pairwise;

        // Temporary terms needed for integer/real demotion/promotion
        real n_numerator_real;
        int nmt_denominator_int;
        real nmt_denominator_real;
        real n_pairwise_real;

        // falling_factorial === (n)_n === N!
        // But falling_factorial requires real values (not ints) so we promote
        // cast our N (an integer) to i (a real)
        n_numerator_real = N;

        // Note that falling_factorial will fail if either value <= 0, so we must
        // max the denominator term to 1. Then promote int j_int to real j_real
        nmt_denominator_int = max(1, N - 2);
        nmt_denominator_real = nmt_denominator_int;

        // We use the temporary term k (real) because falling_factorial returns
        // a real value
        n_pairwise_real = falling_factorial(n_numerator_real, n_numerator_real)
            / (2 * falling_factorial(nmt_denominator_real, nmt_denominator_real));

        // Now convert n_pairwise_real (real) to N_pairwise (int) by looping
        N_pairwise = 1;
        while ((n_pairwise_real + 0.01) > (N_pairwise + 1))
            N_pairwise = N_pairwise + 1;

        return N_pairwise;
    }
}

data {
    int<lower=1> N_nodes;                           // number of estimators (nodes)
    int<lower=1> N_calibrators;                     // number of calibration objects
    int<lower=1> N_visits[N_calibrators, N_nodes];  // number of visits per calibrator.
    int max_visits;

    // Non-spectroscopic (or high fidelity unbiased) measurements
    vector[N_calibrators] mu_calibrator;
    vector<lower=0>[N_calibrators] sigma_calibrator;

    // Spectroscopic estimates 
    //      the estimates provided by the estimators (nodes)
    matrix[N_calibrators, N_nodes, max_visits] estimates;            
    matrix[N_calibrators, N_nodes, max_visits] var_additive;

    // Inverse variance (SNR**-2) of the spectrum
    vector[N_calibrators, N_nodes, max_visits] ivar_spectrum;     
}

transformed data {
    int K; // Number of correlation coefficients between the estimators (nodes)
    K = num_pairwise(N_nodes);
}

parameters {
    // Uncertainty from each estimator
    //      alpha_sq is alpha**2, where \sigma_rand = alpha/SNR
    vector<lower=1>[N_nodes] alpha_sq;
    //      systematic variance in a given estimator (node)
    vector<lower=0>[N_nodes] var_sys_estimator;

    // Correlation coefficients between different nodes.
    vector<lower=-1,upper=+1>[K] rho_estimators;

    // Estimator biases
    vector[N_nodes] bias;

    // God's word
    vector[N_calibrators] truths;
}

model {

    // For each calibrator, calculate the weighted mean from each node.
    for (i in 1:N_calibrators) {
        int a;
        vector mean_mu[N_nodes];
        cov_matrix[N_nodes] Sigma;

        for (j in 1:N_nodes) {
            int K;
            K = N_visits[i, j];

            if (K > 0) {
                // Calculate weights.
                real weights[K];

                weights = 1.0/(alpha_sq[j] * ivar_spectrum[i, j, 1:K]);
                mean_mu[j] = sum(weights * estimates[i, j, 1:K])/sum(weights);
                Sigma[j, j] = 1.0/sum(weights) + var_sys_estimator[j];
            }
        }

        a = 0;
        for (j in 1:N_nodes) {
            for (k in j + 1:N_nodes) {
                Sigma[j, k] = rho_estimators[a] * sqrt(Sigma[j, j] * Sigma[k, k]);
                Sigma[k, j] = rho_estimators[a] * sqrt(Sigma[j, j] * Sigma[k, k]);
                a = a + 1;
            }
            if (N_visits[i, j] == 0) {
                Sigma[j, :] = rep_vector(0.0, N_nodes)';
                Sigma[:, j] = rep_vector(0.0, N_nodes);
                Sigma[j, j] = 1e50;
            }
        }

        mean_mu ~ normal(rep_vector(truths[i], N_nodes) - bias, Sigma);
    }
}

generated quantities {
    // For convenience.
    //real<lower=0> intrinsic_sigma;
    //vector<lower=0>[N_estimators] estimator_sys_sigma;
    //vector<lower=0>[N_estimators] estimator_rand_sigma;

    //intrinsic_sigma = pow(intrinsic_var, 0.5);
}
