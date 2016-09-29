
// Multiple-parameter ensemble model


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
    int<lower=1> N_parameters;
    int<lower=1> N_estimators;
    int<lower=1> N_calibrators;

    // Non-spectroscopic (or high fidelity unbiased) measurements
    vector[N_calibrators * N_parameters] calibrator_mu;
    vector<lower=0>[N_calibrators * N_parameters] calibrator_sigma;

    // Spectroscopic estimates 
    // (allow for multiple estimates of the same star)
    matrix[N_calibrators, N_estimators * N_parameters] mean_estimate;
    matrix[N_calibrators, N_estimators * N_parameters] sqrt_N_estimates_per_estimator;

    // Additive variance is for when no calibrator measurement are finite 
    matrix<lower=0>[N_calibrators, N_estimators * N_parameters] additive_sigma;

}

transformed data {
    int<lower=1> K;
    int<lower=0> N_rho_parameter_terms;
    int<lower=0> N_rho_estimator_terms;

    matrix<lower=0>[N_estimators * N_parameters, N_estimators * N_parameters] additive_var[N_calibrators];

    K = N_estimators * N_parameters;
    N_rho_parameter_terms = num_pairwise(N_parameters);
    N_rho_estimator_terms = N_parameters * num_pairwise(N_estimators);

    for (i in 1:N_calibrators)
        additive_var[i] = rep_matrix(additive_sigma[i], K) .* rep_matrix(additive_sigma[i], K)';

}


parameters {
    // Intrinsic uncertainty in the model
    //vector<lower=0>[N_parameters] intrinsic_var;

    // Uncertainty from each estimator (node), both random and systematic
    matrix<lower=0>[N_estimators, N_parameters] estimator_sys_sigma;
    //matrix<lower=0>[N_estimators, N_parameters] estimator_rand_sigma;

    // God's word
    matrix[N_calibrators, N_parameters] truths;

    // Correlations between parameters.
    vector<lower=-1, upper=1>[N_rho_parameter_terms] rho_parameters;
    vector<lower=-1, upper=1>[N_rho_estimator_terms] rho_estimators;
}

transformed parameters {

    // Build our fuck-off-big covariance matrix

    cov_matrix[K] Sigma[N_calibrators];
    {
        matrix[K, K] rho;
        matrix[K, K] estimator_sys_variance;

        // Initialize rho as an identity matrix
        rho = diag_matrix(rep_vector(1.0, K));
        {
            // Temporary local variables to make code cleaner
            int a;
            int b;
            int offset;
            matrix[N_estimators, N_estimators] rep_rho_parameters;
            a = 1;
            b = 1;

            for (i in 1:N_parameters) {
                // rho_{estimator-a,estimator-b} terms
                for (j in 1:N_estimators) {
                    for (k in j + 1:N_estimators) {
                        offset = (i - 1) * N_estimators;

                        rho[j + offset, k + offset] = rho_estimators[a];
                        rho[k + offset, j + offset] = rho_estimators[a];
                        a = a + 1;
                    }
                }

                // rho_{parameter-a,parameter-b} terms
                for (j in i + 1:N_parameters) {

                    rep_rho_parameters = rep_matrix(
                        rho_parameters[b], N_estimators, N_estimators);

                    rho[1 + N_estimators * (i - 1):N_estimators * i,
                        1 + N_estimators * (j - 1):N_estimators * j] = rep_rho_parameters;
                    rho[1 + N_estimators * (j - 1):N_estimators * j,
                        1 + N_estimators * (i - 1):N_estimators * i] = rep_rho_parameters;
                    b = b + 1;
                }
            }

        }
        
        estimator_sys_variance = 
            (rep_matrix(to_vector(estimator_sys_sigma), K)
                .* rep_matrix(to_vector(estimator_sys_sigma), K)');
    
        for (i in 1:N_calibrators) {
            Sigma[i] = rho .* (
                // Systematic uncertainty from a single node
                
                estimator_sys_variance
                
                + 
                
                // Random uncertainty using multiple estimates from a single node 
                /*
                (
                    (rep_matrix(to_vector(estimator_rand_sigma) ./ to_vector(sqrt_N_estimates_per_estimator[i]), K) 
                        .* rep_matrix(to_vector(estimator_rand_sigma) ./ to_vector(sqrt_N_estimates_per_estimator[i]), K)')
                )
                +
                */

                // Additive variance due to missing data values
                + additive_var[i]
                
            );  
        }     
    }

}

model {
    calibrator_mu ~ normal(to_vector(truths), calibrator_sigma); 
    for (i in 1:N_calibrators)
        mean_estimate[i] ~ multi_normal(
            to_vector(rep_matrix(truths[i], N_estimators)),
            Sigma[i]);
}

generated quantities {
    // For convenience.
    //real<lower=0> intrinsic_sigma;
    //vector<lower=0>[N_estimators] estimator_sys_sigma;
    //vector<lower=0>[N_estimators] estimator_rand_sigma;

    //intrinsic_sigma = pow(intrinsic_var, 0.5);
}