
// Single-parameter ensemble model with correlation coefficients between nodes

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
    int<lower=1> N_estimators;          // number of estimators (nodes)
    int<lower=1> N_calibrators;         // number of calibration objects
    int<lower=1> N_calibrator_visits;   // number of measurements (or visits) of calibration objects

    // Non-spectroscopic (or high fidelity unbiased) measurements
    //      the calibrator parameters (e.g., teff determined by non-spectroscopic method)
    vector[N_calibrators] mu_calibrator;
    //      the uncertainty in the calibrator values
    vector<lower=0>[N_calibrators] sigma_calibrator;

    // Spectroscopic estimates 
    //      the estimates provided by the estimators (nodes)
    matrix[N_calibrator_visits, N_estimators] estimates;            
    //      the index of the calibrator corresponding to each visit
    int calibrator_index[N_calibrator_visits];                      

    // Additive variance is for when no calibrator measurement are finite 
    vector<lower=0>[N_estimators] var_additive[N_calibrator_visits]; 

    // Inverse variance (SNR**-2) of the spectrum
    vector[N_calibrator_visits] ivar_spectrum;     

    matrix[N_estimators, N_estimators] Sigma_c0[N_calibrator_visits];
    matrix[N_estimators, N_estimators] Sigma_c1[N_calibrator_visits];                     
}

transformed data {
    int K; // Number of correlation coefficients between the estimators (nodes)
    K = num_pairwise(N_estimators);
}

parameters {
    // Intrinsic uncertainty in the model
    real<lower=0> var_intrinsic;

    // Uncertainty from each estimator
    //      alpha_sq is alpha**2, where \sigma_rand = alpha/SNR
    vector<lower=100>[N_estimators] alpha_sq;
    //      systematic variance in a given estimator (node)
    vector<lower=0>[N_estimators] var_sys_estimator;

    // Correlation coefficients between different nodes.
    vector<lower=-1,upper=+1>[K] rho_estimators;

    // Estimator biases
    vector[N_estimators] c0_estimators;

    // God's word
    vector[N_calibrators] truths;
}

transformed parameters {
    // This is where the pain begins.
    cov_matrix[N_estimators] Sigma[N_calibrator_visits];
    {
        int a;
        real sigma_j;
        real sigma_k;
        for (i in 1:N_calibrator_visits) {
            a = 1;
            for (j in 1:N_estimators) {
                for (k in j:N_estimators) {
                    sigma_j = sqrt(
                        var_intrinsic + 
                        var_sys_estimator[j] + 
                        alpha_sq[j] * ivar_spectrum[i]);

                    sigma_k = sqrt(
                        var_intrinsic +
                        var_sys_estimator[k] + 
                        alpha_sq[k] * ivar_spectrum[i]);

                    if (j == k) {
                        Sigma[i,j,j] = sigma_j * sigma_j;
                    }
                    else {
                        Sigma[i,j,k] = rho_estimators[a] * sigma_j * sigma_k;
                        Sigma[i,k,j] = rho_estimators[a] * sigma_j * sigma_k; 
                        a = a + 1;
                    }
                }
            }

            for (j in 1:N_estimators) {
                if (var_additive[i, j] > 1) {
                    Sigma[i, j, :] = rep_vector(0.0, N_estimators)';
                    Sigma[i, :, j] = rep_vector(0.0, N_estimators);
                    Sigma[i, j, j] = var_additive[i, j];
                }
            }

        }
    }
}

model {
    mu_calibrator ~ normal(to_vector(truths), sigma_calibrator); 
    for (i in 1:N_calibrator_visits) {
        estimates[i] ~ multi_normal(
            to_vector(rep_vector(truths[calibrator_index[i]], N_estimators)) - c0_estimators,
            Sigma[i]);
    }
}

generated quantities {
    // For convenience.
    //real<lower=0> intrinsic_sigma;
    //vector<lower=0>[N_estimators] estimator_sys_sigma;
    //vector<lower=0>[N_estimators] estimator_rand_sigma;

    //intrinsic_sigma = pow(intrinsic_var, 0.5);
}
