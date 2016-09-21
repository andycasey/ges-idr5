
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

    // Non-spectroscopic measurements
    matrix[N_calibrators, N_parameters] calibrator_mu;
    matrix<lower=0>[N_calibrators, N_parameters] calibrator_sigma;

    // Spectroscopic measurements (allow for multiple measurements of the same
    // star)
    matrix[N_estimators, N_parameters] mean_estimate[N_calibrators];
    matrix[N_estimators, N_parameters] N_estimates_per_estimator[N_calibrators];

    // Additive variance is for when no calibrator measurement are finite 
    vector<lower=0>[N_estimators, N_parameters] additive_var[N_calibrators];
}

transformed data {
    int<lower=0> N_rho_parameter_terms;
    int<lower=0> N_rho_estimator_terms;

    N_rho_parameter_terms = num_pairwise(N_parameters);
    N_rho_estimator_terms = num_pairwise(N_estimators);
}

parameters {
    // Intrinsic uncertainty in the model
    //real<lower=0> intrinsic_var;

    // Uncertainty from each estimator (node), both random and systematic
    matrix<lower=0>[N_estimators, N_parameters] estimator_sys_var;
    matrix<lower=0>[N_estimators, N_parameters] estimator_rand_var;

    // God's word
    matrix[N_calibrators, N_parameters] truths;

    // Correlations between parameters.
    vector<lower=-1, upper=1>[N_rho_parameter_terms] rho_parameters;
    vector<lower=-1, upper=1>[N_rho_estimator_terms] rho_estimators;
}

transformed parameters {
    // Build our fuck-off big covariance matrix
    matrix[N_estimators, N_estimators] covariance[N_calibrators];
    for (i in 1:N_calibrators)
        covariance[i] <- rep_matrix(intrinsic_var, N_estimators, N_estimators)
             + diag_matrix(estimator_sys_var)
             + diag_matrix(estimator_rand_var ./ N_estimates_per_estimator[i])
             + diag_matrix(additive_var[i]);
}

model {
    calibrator_mu ~ normal(truths, calibrator_sigma); 
    for (i in 1:N_calibrators) 
        mean_estimate[i] ~ multi_normal(
            rep_vector(truths[i], N_estimators), covariance[i]);
}
