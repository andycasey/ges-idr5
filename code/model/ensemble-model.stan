
// Simple ensemble model with multiple measurements

data {
    int<lower=1> N_estimators;
    int<lower=1> N_calibrators;
    // TEST
    int<lower=1> N_parameters;

    // Non-spectroscopic measurements
    vector[N_calibrators] calibrator_mu;
    vector<lower=0>[N_calibrators] calibrator_sigma;

    // Spectroscopic measurements (allow for multiple measurements of the same
    // star)
    vector[N_estimators] mean_estimate[N_calibrators];
    vector[N_estimators] N_estimates_per_estimator[N_calibrators];

    // Additive variance is for when no calibrator measurement are finite 
    vector<lower=0>[N_estimators] additive_var[N_calibrators];
}

parameters {
    // Intrinsic uncertainty in the model
    real<lower=0> intrinsic_var;

    // Uncertainty from each estimator (node), both random and systematic
    vector<lower=0>[N_estimators] estimator_sys_var;
    vector<lower=0>[N_estimators] estimator_rand_var;

    // God's word
    vector[N_calibrators] truths;

    // Test 
    // Number of pair-wise correlation coeficitions is N!/(2(N-2)!) where N is
    // the number of parameters and N! is the falling factorial such that
    // N! = (N)_N
}

transformed parameters {
    matrix[N_estimators, N_estimators] covariance[N_calibrators];
    for (i in 1:N_calibrators)
        covariance[i] = rep_matrix(intrinsic_var, N_estimators, N_estimators)
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
/*
generated quantities {
    // For convenience.
    real<lower=0> intrinsic_sigma;
    vector<lower=0>[N_estimators] estimator_sys_sigma;
    vector<lower=0>[N_estimators] estimator_rand_sigma;

    intrinsic_sigma = pow(intrinsic_var, 0.5);
    estimator_sys_sigma <- pow(estimator_sys_var, 0.5);
    estimator_rand_sigma <- pow(estimator_rand_var, 0.5);
}
*/