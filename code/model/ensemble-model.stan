
// Simple ensemble model with multiple measurements
// Dreamt up by Sergey Koposov, shittily extended upon by Andy Casey

data {
    int<lower=1> N_nodes;
    int<lower=1> N_benchmarks;

    // Non-spectroscopic measurements
    vector[N_benchmarks] non_spec_measurement;
    vector<lower=0>[N_benchmarks] non_spec_sigma;

    // Spectroscopic measurements (allow for multiple measurements of the same
    // star)
    vector[N_nodes] averaged_node_measurement[N_benchmarks];
    vector[N_nodes] number_of_node_measurements[N_benchmarks];

    // Additive variance is for when no benchmark measurement are finite 
    vector<lower=0>[N_nodes] averaged_additive_variance[N_benchmarks];
}

parameters {
    // Uncertainty in the data
    real<lower=0> var_intrinsic;

    // Uncertainty in each node
    vector<lower=0>[N_nodes] var_node_sys;
    vector<lower=0>[N_nodes] var_node_rand;

    // God's word
    vector[N_benchmarks] truths;
}

transformed parameters {
    matrix[N_nodes, N_nodes] covariance[N_benchmarks];
    for (i in 1:N_benchmarks)
        covariance[i] <- rep_matrix(var_intrinsic, N_nodes, N_nodes)
             + diag_matrix(var_node_sys)
             + diag_matrix(var_node_rand ./ number_of_node_measurements[i]) # element-wise division
             + diag_matrix(averaged_additive_variance[i]);
}

model {
    non_spec_measurement ~ normal(truths, non_spec_sigma); 
    for (i in 1:N_benchmarks) 
        averaged_node_measurement[i] ~ multi_normal(
            rep_vector(truths[i], N_nodes), covariance[i]);
}
