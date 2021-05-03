data { 
    // Dimensional parameters
    int<lower=1> J; // Number of biological replicates for yield
    int<lower=1> N_yield; // Number of yield measurements
    int<lower=1> N_calib; // Number of calibration measurements
    int<lower=1, upper=J> idx[N_yield]; // ID vector for replicates

    // Calibration observables
    vector<lower=0>[N_calib] calib_conc; 
    vector<lower=0>[N_calib] calib_rel_areas;

    // Yield observables
    vector<lower=0>[N_yield] optical_density;
    vector<lower=0>[N_yield] yield_rel_areas; 
}

transformed data {
    vector[N_calib] centered_conc = (calib_conc - mean(calib_conc)) / sd(calib_conc);
    vector[N_calib] centered_rel_area = (calib_rel_areas - mean(calib_rel_areas)) / sd(calib_rel_areas);
}

parameters { 
    // Calibration parameters
    real calib_slope_cent;
    real calib_inter_cent;
    real<lower=0> calib_sigma_cent;

    // Hyperparameters
    real yield_inter_area_mu;
    real<lower=0> yield_inter_area_sigma;
    real yield_slope_area_mu;
    real<lower=0> yield_slope_area_sigma;
    
    // Low level parameters
    vector[J] yield_inter_area;
    vector[J] yield_slope_area;

    // Homoscedastic error
    real<lower=0> yield_area_sigma; 
}

transformed parameters { 
    real calib_slope = calib_slope_cent * sd(calib_rel_areas)/sd(calib_conc);
    real calib_inter = sd(calib_rel_areas) * (calib_inter_cent - calib_slope_cent * mean(calib_conc) / sd(calib_conc)) + mean(calib_rel_areas);
    real calib_sigma = sd(calib_rel_areas) * calib_sigma_cent;
}

model {

    // Calibration priors
    calib_inter_cent ~ std_normal();
    calib_slope_cent ~ std_normal();
    calib_sigma_cent ~ std_normal();

    // Hyper parameters
    yield_inter_area_mu ~ std_normal();
    yield_slope_area_mu ~ std_normal();
    yield_inter_area_sigma ~ normal(0, 0.1);
    yield_slope_area_sigma ~ normal(0, 0.1);

    // Low-level priors
    yield_inter_area ~ normal(yield_inter_area_mu, yield_inter_area_sigma);
    yield_slope_area ~ normal(yield_slope_area_mu, yield_slope_area_sigma);

    // Homoscedastic error
    yield_area_sigma ~ std_normal();

    // Calibration likelihood
    centered_rel_area ~ normal(calib_slope_cent * centered_conc + calib_inter_cent, calib_sigma_cent);

    // Yield likelihood
    yield_rel_areas ~ normal(yield_slope_area[idx] .* optical_density + yield_inter_area[idx], yield_area_sigma);
}

generated quantities {
    // Compute the concentration for easy accesssion    
    vector[N_yield] yield_concs = (yield_rel_areas - calib_inter) / calib_slope;
    vector[J] yield_inter = (yield_inter_area - calib_inter) / calib_slope;
    vector[J] yield_slope = (yield_slope_area - calib_inter) / calib_slope;
    real yield_inter_mu = (yield_inter_area_mu - calib_inter) / calib_slope;
    real yield_slope_mu = (yield_slope_area_mu - calib_inter) / calib_slope;

}