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

parameters { 
    // Calibration parameters
    real<lower=0> calib_slope;
    real<lower=0> calib_inter;
    real<lower=0> calib_sigma;

    // Hyperparameters
    real<lower=0> yield_inter_mu;
    real<lower=0> yield_inter_sigma;
    real yield_slope_mu;
    real<lower=0> yield_slope_sigma;
    
    // Low level parameters
    vector[J] yield_inter;
    vector[J] yield_slope;

    // Homoscedastic error
    real<lower=0> yield_sigma; 
}


model {
    // Calculate the concentrations given the calibration inference
    vector[N_yield] yield_concs = (yield_rel_areas - calib_inter) / calib_slope;

    // Calibration priors
    calib_inter ~ std_normal();
    calib_slope ~ std_normal();
    calib_sigma ~ std_normal();

    // Hyper parameters
    yield_inter_mu ~ std_normal();
    yield_slope_mu ~ std_normal();
    yield_inter_sigma ~ normal(0, 0.1);
    yield_slope_sigma ~ normal(0, 0.1);

    // Low-level priors
    yield_inter ~ normal(yield_inter_mu, yield_inter_sigma);
    yield_slope ~ normal(yield_slope_mu, yield_slope_sigma);

    // Homoscedastic error
    yield_sigma ~ std_normal();

    // Calibration likelihood
    calib_rel_areas ~ normal(calib_slope * calib_conc + calib_inter, calib_sigma);

    // Yield likelihood
    yield_concs ~ normal(yield_slope[idx] .* optical_density + yield_inter[idx], yield_sigma);
}

generated quantities {
    // Compute the concentration for easy accesssion    
    vector[N_yield] yield_concs = (yield_rel_areas - calib_inter) / calib_slope;
}