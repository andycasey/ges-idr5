
// Single-parameter ensemble model with correlation coefficients between nodes

data {
    int<lower=1> N; // number of nodes (estimators)
    int<lower=1> C; // number of calibrators (benchmarks)
    int<lower=1> V; // maximum number of visits to any calibrator
    int<lower=1> visits[C]; // number of visits to each calibrator
    int<lower=1> S; // the number of astrophysical parameters that can contribute
                    // to the node systematic variance (rank-ordered per param)
    
    int<lower=0> TM; // total number of missing data points.
    int is_missing[C, N, V]; // whether a node estimate is missing (>0 = missing)

    matrix[N, V] estimates[C]; // estimates of the stellar property
    matrix[V, C] spectrum_isnr; // inverse-snr of the spectrum (1/SNR)  

    real mu_calibrator[C]; // the non-spectroscopic or high fidelity calibrator value
    real sigma_calibrator[C]; // the 1\sigma uncertainty on the calibrator value

    real all_mu_calibrator[C, S]; // all non-spectroscopic calibrator values 
                                  // available (for modeling systematic variance)

    real lower_alpha_sq; // lower bound on alpha parameter for f(snr)
    real upper_alpha_sq; // upper bound on alpha parameter for f(snr)

    real lower_bound; // lower univariate bound on the missing data values
    real upper_bound; // upper univariate bound on the missing data values
}

transformed data {
    int TCV; // total number of calibrator visits
    matrix[9, C] AMCS; // a design matrix of all_mu_calibrator values scaled
                       // between (0, 1)
                       // #MAGIC HACK: 9 is the number of terms in a quadratic
                       //              3-parameter system with a separated constant

    TCV = sum(visits);
    for (c in 1:C) {
        real teff_scale;
        real teff_offset;
        real teff_norm;
        
        real logg_scale;
        real logg_offset;
        real logg_norm;

        real feh_scale;
        real feh_offset;
        real feh_norm;

        teff_offset = min(all_mu_calibrator[:, 1]);
        teff_scale = max(all_mu_calibrator[:, 1]) - teff_offset;
        teff_norm = 1.0 - (all_mu_calibrator[c, 1] - teff_offset)/teff_scale;

        logg_offset = min(all_mu_calibrator[:, 2]);
        logg_scale = max(all_mu_calibrator[:, 2]) - logg_offset;
        logg_norm = 1.0 - (all_mu_calibrator[c, 2] - logg_offset)/logg_scale;
        
        feh_offset = min(all_mu_calibrator[:, 3]);
        feh_scale = max(all_mu_calibrator[:, 3]) - feh_offset;
        feh_norm = 1.0 - (all_mu_calibrator[c, 3] - feh_offset)/feh_scale;

        AMCS[1, c] = pow(teff_norm, 2);
        AMCS[2, c] = pow(logg_norm, 2);
        AMCS[3, c] = pow(feh_norm, 2);

        AMCS[4, c] = teff_norm;
        AMCS[5, c] = logg_norm;
        AMCS[6, c] = feh_norm;

        AMCS[7, c] = teff_norm * logg_norm;
        AMCS[8, c] = teff_norm * feh_norm;
        AMCS[9, c] = logg_norm * feh_norm;
    }
}

parameters {
    real truths[C]; // God's word: true values of the calibrators.
    real biases[N]; // biases in the individual nodes
    real<lower=lower_bound, upper=upper_bound> missing_estimates[TM];

    // Cholesky factor of a correlation matrix
    cholesky_factor_corr[N] L_corr;
    
    // \alpha:  parameter to model uncertainty as a function of SNR s.t.
    //          \_variance_{random} = \alpha^2 * spectrum_isnr
    vector<lower=lower_alpha_sq, upper=upper_alpha_sq>[N] alpha_sq; 
    vector<lower=0, upper=upper_alpha_sq>[N] vs_c; // constant systematic term
    //matrix[9, N] vs_theta; // coefficients to increase relative systematic
    //                       // uncertainty across parameter space.
    //                       // #MAGIC: 9 matches AMCS for the same reason.

    real vs_tb1;
    real vs_tb2;
    real vs_tb3;
    real vs_tb4;
    real vs_tb5;
    real vs_tb6;
    real vs_tb7;
    real vs_tb8;

    real vs_lb1;
    real vs_lb2;
    real vs_lb3;
    real vs_lb4;
    real vs_lb5;
    real vs_lb6;
    real vs_lb7;
    real vs_lb8;

    real vs_fb1;
    real vs_fb2;
    real vs_fb3;
    real vs_fb4;
    real vs_fb5;
    real vs_fb6;
    real vs_fb7;
    real vs_fb8;

    real<lower=pow(vs_tb1/2.0, 2)> vs_ta1;
    real<lower=pow(vs_tb2/2.0, 2)> vs_ta2;
    real<lower=pow(vs_tb3/2.0, 2)> vs_ta3;
    real<lower=pow(vs_tb4/2.0, 2)> vs_ta4;
    real<lower=pow(vs_tb5/2.0, 2)> vs_ta5;
    real<lower=pow(vs_tb6/2.0, 2)> vs_ta6;
    real<lower=pow(vs_tb7/2.0, 2)> vs_ta7;
    real<lower=pow(vs_tb8/2.0, 2)> vs_ta8;

    real<lower=pow(vs_lb1/2.0, 2)> vs_la1;
    real<lower=pow(vs_lb2/2.0, 2)> vs_la2;
    real<lower=pow(vs_lb3/2.0, 2)> vs_la3;
    real<lower=pow(vs_lb4/2.0, 2)> vs_la4;
    real<lower=pow(vs_lb5/2.0, 2)> vs_la5;
    real<lower=pow(vs_lb6/2.0, 2)> vs_la6;
    real<lower=pow(vs_lb7/2.0, 2)> vs_la7;
    real<lower=pow(vs_lb8/2.0, 2)> vs_la8;

    real<lower=pow(vs_fb1/2.0, 2)> vs_fa1;
    real<lower=pow(vs_fb2/2.0, 2)> vs_fa2;
    real<lower=pow(vs_fb3/2.0, 2)> vs_fa3;
    real<lower=pow(vs_fb4/2.0, 2)> vs_fa4;
    real<lower=pow(vs_fb5/2.0, 2)> vs_fa5;
    real<lower=pow(vs_fb6/2.0, 2)> vs_fa6;
    real<lower=pow(vs_fb7/2.0, 2)> vs_fa7;
    real<lower=pow(vs_fb8/2.0, 2)> vs_fa8;

    real<lower=0> vs_tc7[N]; // Cross term coefficients for teff * logg 
    real<lower=0> vs_tc8[N]; // Cross term coefficients for teff * feh
    real<lower=0> vs_tc9[N]; // Cross term coefficients for logg * feh 

}

transformed parameters {
    cov_matrix[N] Sigma[TCV]; // covariance matrix
    matrix[N, V] full_rank_estimates[C]; // array containing known (data) and 
                                         // unknown (parameter) estimates
    vector<lower=-3>[N] relative_sigma_sys; // relative increase in systematic
                                            // uncertainty for a given node
    {
        int mi; // an index for unknown data points
        int vi; // an index for total visits
        mi = 1;
        vi = 1;

        for (c in 1:C) { // For each calibrator

            // Construct the relative systematic uncertainty here
            relative_sigma_sys[1] = vs_ta1 * AMCS[1, c]
                                  + vs_la1 * AMCS[2, c]
                                  + vs_fa1 * AMCS[3, c]

                                  + vs_tb1 * AMCS[4, c]
                                  + vs_lb1 * AMCS[5, c]
                                  + vs_fb1 * AMCS[6, c]

                                  + vs_tc7[1] * AMCS[7, c]
                                  + vs_tc8[1] * AMCS[8, c]
                                  + vs_tc9[1] * AMCS[9, c];

            relative_sigma_sys[2] = vs_ta2 * AMCS[1, c]
                                  + vs_la2 * AMCS[2, c]
                                  + vs_fa2 * AMCS[3, c]

                                  + vs_tb2 * AMCS[4, c]
                                  + vs_lb2 * AMCS[5, c]
                                  + vs_fb2 * AMCS[6, c]

                                  + vs_tc7[2] * AMCS[7, c]
                                  + vs_tc8[2] * AMCS[8, c]
                                  + vs_tc9[2] * AMCS[9, c];

            relative_sigma_sys[3] = vs_ta3 * AMCS[1, c]
                                  + vs_la3 * AMCS[2, c]
                                  + vs_fa3 * AMCS[3, c]

                                  + vs_tb3 * AMCS[4, c]
                                  + vs_lb3 * AMCS[5, c]
                                  + vs_fb3 * AMCS[6, c]

                                  + vs_tc7[3] * AMCS[7, c]
                                  + vs_tc8[3] * AMCS[8, c]
                                  + vs_tc9[3] * AMCS[9, c];                                  

            relative_sigma_sys[4] = vs_ta4 * AMCS[1, c]
                                  + vs_la4 * AMCS[2, c]
                                  + vs_fa4 * AMCS[3, c]

                                  + vs_tb4 * AMCS[4, c]
                                  + vs_lb4 * AMCS[5, c]
                                  + vs_fb4 * AMCS[6, c]

                                  + vs_tc7[4] * AMCS[7, c]
                                  + vs_tc8[4] * AMCS[8, c]
                                  + vs_tc9[4] * AMCS[9, c];

            relative_sigma_sys[5] = vs_ta5 * AMCS[1, c]
                                  + vs_la5 * AMCS[2, c]
                                  + vs_fa5 * AMCS[3, c]

                                  + vs_tb5 * AMCS[4, c]
                                  + vs_lb5 * AMCS[5, c]
                                  + vs_fb5 * AMCS[6, c]

                                  + vs_tc7[5] * AMCS[7, c]
                                  + vs_tc8[5] * AMCS[8, c]
                                  + vs_tc9[5] * AMCS[9, c];

            relative_sigma_sys[6] = vs_ta6 * AMCS[1, c]
                                  + vs_la6 * AMCS[2, c]
                                  + vs_fa6 * AMCS[3, c]

                                  + vs_tb6 * AMCS[4, c]
                                  + vs_lb6 * AMCS[5, c]
                                  + vs_fb6 * AMCS[6, c]

                                  + vs_tc7[6] * AMCS[7, c]
                                  + vs_tc8[6] * AMCS[8, c]
                                  + vs_tc9[6] * AMCS[9, c];

            relative_sigma_sys[7] = vs_ta7 * AMCS[1, c]
                                  + vs_la7 * AMCS[2, c]
                                  + vs_fa7 * AMCS[3, c]

                                  + vs_tb7 * AMCS[4, c]
                                  + vs_lb7 * AMCS[5, c]
                                  + vs_fb7 * AMCS[6, c]

                                  + vs_tc7[7] * AMCS[7, c]
                                  + vs_tc8[7] * AMCS[8, c]
                                  + vs_tc9[7] * AMCS[9, c];

            relative_sigma_sys[8] = vs_ta8 * AMCS[1, c]
                                  + vs_la8 * AMCS[2, c]
                                  + vs_fa8 * AMCS[3, c]

                                  + vs_tb8 * AMCS[4, c]
                                  + vs_lb8 * AMCS[5, c]
                                  + vs_fb8 * AMCS[6, c]

                                  + vs_tc7[8] * AMCS[7, c]
                                  + vs_tc8[8] * AMCS[8, c]
                                  + vs_tc9[8] * AMCS[9, c];


            //for (n in 1:N)
            //    relative_sigma_sys[n, c] = sum(AMCS[:, c] .* vs_theta[:, n]);

            // For each visit of the calibrator
            for (v in 1:visits[c]) {
                vector[N] sigma;

                // Build the diagonal of the covariance matrix
                for (n in 1:N) {

                    // The total (random and systematic) uncertainty
                    sigma[n] = sqrt(
                        alpha_sq[n] * spectrum_isnr[v, c] +
                        vs_c[n] * (4.0 + relative_sigma_sys[n]));

                    if (is_missing[c, n, v]) {
                        full_rank_estimates[c, n, v] = missing_estimates[mi] - biases[n];
                        mi = mi + 1;
                    }
                    else {
                        full_rank_estimates[c, n, v] = estimates[c, n, v] - biases[n];
                    }
                }
                // Construct the covariance matrix from Cholesky factors
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
