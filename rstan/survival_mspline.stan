functions {
  // -- Which entries are equal to a given value --
  array[] int which(array[] int x, int y) {
    int n = num_elements(x);
    array[n] int w; // Over-allocate w and then truncate later
    int c = 1;
    for (i in 1:n) {
      if (x[i] == y) {
        w[c] = i;
        c += 1;
      }
    }
    return w[1:(c-1)];
  } 
    
  array[] int nwhich_all(array[] int x, int max_id) {
    array[max_id] int w = rep_array(0, max_id);
    for (i in 1:num_elements(x)) {
      if (x[i]) w[x[i]] += 1;
    }
    return w;
  }
  
  // -- M-splines log-Survival and log-Hazard --
  // IPD: Evaluate for all event times within each study
  vector lS (matrix ibasis, vector eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = vector of spline coefficients
    // eta = log rate (linear predictor)
    return - (ibasis * scoef) .* exp(eta);
  }

  vector lh (matrix basis, vector eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(basis * scoef) + eta;
  }
  
  // AgD: Evaluate for a vector of integration points (theta/eta) for a single event time
  vector lS_a (row_vector ibasis, vector eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = vector of spline coefficients
    // eta = log rate (linear predictor)
    return - (ibasis * scoef) * exp(eta);
  }

  vector lh_a (row_vector basis, vector eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(basis * scoef) + eta;
  }
  
  // -- Log likelihood functions --
  // IPD version
  vector loglik(matrix time,    // Basis evaluated at event/cens time
                matrix itime,   // Integrated basis evaluated at event/cens time
                array[] int status,   // Censoring status
                vector eta,     // Linear predictor
                vector scoef) { // Spline coefficients

    vector[num_elements(eta)] l;
    int nw1 = sum(status);
    array[nw1] int w1 = which(status, 1);

    // All observations get survival contribution
    l = lS(itime, eta, scoef);

    // Observed events get additional hazard contribution
    if (nw1) {
      l[w1] += lh(time[w1], eta[w1], scoef);
    }

    return l;
  }
  
  // AgD version
  vector loglik_a(row_vector time,  // Basis evaluated at event/cens time
                  row_vector itime, // Integrated basis evaluated at event/cens time
                  int status,       // Censoring status
                  vector eta,       // Linear predictor
                  vector scoef) {   // Spline coefficients
    vector[num_elements(eta)] l;
  
    if (status == 0) { // Right censored
      l = lS_a(itime, eta, scoef);
    } else if (status == 1) { // Observed
      l = lS_a(itime, eta, scoef) + lh_a(time, eta, scoef);
    } 
  
    return l;
  }
}
data {
  // -- Constants --
  int<lower=0> ns_ipd; // number of IPD studies
  int<lower=0> ns_agd; // number of AgD studies
  int<lower=0> ni_ipd; // total number of IPD individuals
  int<lower=0> ni_agd; // total number of AgD data points
  int<lower=1> nt; // number of treatments
  int<lower=1> nint; // number of samples for numerical integration
  int<lower=0> nX; // number of columns in design matrix
  
  // -- Survival data --
  array[ni_ipd] int<lower=0, upper=1> ipd_status; // event status (0 = censored, 1 = event)
  array[ni_agd] int<lower=0, upper=1> agd_status;
  
  int<lower=2> n_scoef; // Number of spline coefficients
  matrix[ni_ipd, n_scoef] ipd_time;   // M-spline basis evaluated at event/censoring times
  matrix[ni_agd, n_scoef] agd_time;
  matrix[ni_ipd, n_scoef] ipd_itime;  // Integrated M-spline basis evaluated at event/censoring times
  matrix[ni_agd, n_scoef] agd_itime;
  
  array[ni_ipd + ni_agd] int<lower=1> aux_id; // aux id (e.g. study id for a PH/AFT model)
  
  // -- Integration point details --
  array[ni_agd] int<lower=1> int_id; // integration sample id (i.e. for an arm or study)
  
  // -- Design matrix --
  // Either X, or QR Q if QR=TRUE
  matrix[ni_ipd + (ni_agd ? nint * max(int_id) : 0), nX] X;
  
  // -- Thin QR decomposition --
  int<lower=0, upper=1> QR; // QR decomposition flag (1=used QR decomposition)
  matrix[nX, nX] R_inv;
  
  // -- Priors --
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_reg_sd;
  real<lower=0> prior_trt_sd;
  array[max(aux_id)] vector[n_scoef-1] prior_lscoef_location;  // prior logistic mean
  array[max(aux_id)] vector[n_scoef-1] prior_lscoef_weight;  // rw1 prior weights
  real<lower=0> prior_sigma_sd; // rw1 smoothing sd
}
transformed data {
  // Total number of studies
  int ns = ns_ipd + ns_agd;
  
  // Number of aux coefficient vectors
  int n_aux = max(aux_id);
  
  // Split aux_id by IPD and AgD rows
  array[ni_ipd] int<lower=1> aux_id_ipd = aux_id[1:ni_ipd];
  array[ni_agd] int<lower=1> aux_id_agd = aux_id[(ni_ipd + 1):(ni_ipd + ni_agd)];

  // Split design matrix into IPD and AgD
  matrix[ni_ipd, nX] X_ipd;
  matrix[ni_agd ? nint * max(int_id) : 0, nX] X_agd;

  // For efficiency we will loop over aux_id's, vectorising within each group
  // Number of observations within each ID
  array[n_aux] int ni_aux_id_ipd = nwhich_all(aux_id_ipd, n_aux);
  array[n_aux] int ni_aux_id_agd = nwhich_all(aux_id_agd, n_aux);
  
  // Index lookup - which observations belong to group i?
  array[n_aux, max(ni_aux_id_ipd)] int wi_aux_id_ipd;
  array[n_aux, max(ni_aux_id_agd)] int wi_aux_id_agd;

  for (i in 1:n_aux) {
    if (ni_aux_id_ipd[i]) wi_aux_id_ipd[i, 1:ni_aux_id_ipd[i]] = which(aux_id_ipd, i);
    if (ni_aux_id_agd[i]) wi_aux_id_agd[i, 1:ni_aux_id_agd[i]] = which(aux_id_agd, i);
  }
  
  if (ni_ipd) X_ipd = X[1:ni_ipd, ];
  if (ni_agd) X_agd = X[(ni_ipd + 1):(ni_ipd + nint * max(int_id)), ];
  
}
parameters {
  // Parameters on QR scale
  vector[nX] beta_tilde;
  
  // Spline shrinkage SDs and non-centered parameterisation smoothing
  vector<lower=0>[n_aux] sigma;
  array[n_aux] vector[n_scoef-1] u_aux;
}
transformed parameters {
  // -- log likelihood --
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd] log_L_agd;
  
  // -- Spline coefficients --
  array[n_aux] vector[n_scoef-1] lscoef; // logit spline coefficients
  array[n_aux] vector[n_scoef] scoef; // spline coefficients

  // -- Back-transformed parameters --
  vector[nX] allbeta = QR? R_inv * beta_tilde : beta_tilde;
  // Study baselines
  vector[ns] mu = allbeta[1:ns];
  // Treatment effects
  vector[nt - 1] gamma = allbeta[(ns + 1):(ns + nt - 1)];
  // Regression terms
  vector[nX - ns - (nt - 1)] beta = allbeta[(ns + nt):nX];

  // -- Construct spline coefficients with random walk prior around constant hazard --
  for (i in 1:n_aux) {
    lscoef[i] = cumulative_sum(u_aux[i] .* prior_lscoef_weight[i]) * sigma[i] + prior_lscoef_location[i];
    scoef[i] = softmax(append_row(0, lscoef[i]));
  }

  // -- IPD log likelihood --
  // Loop over aux parameters (i.e. by study)
  if (ni_ipd) {
    vector[ni_ipd] eta_ipd = X_ipd * beta_tilde;
    
    for (i in 1:n_aux) {
      int ni = ni_aux_id_ipd[i];
  
      if (ni) {
        array[ni] int wi = wi_aux_id_ipd[i, 1:ni];
        log_L_ipd[wi] = loglik(ipd_time[wi],
                               ipd_itime[wi],
                               ipd_status[wi],
                               eta_ipd[wi],
                               scoef[i]);
      }
    }
  }

  // -- AgD log likelihood --
  if (ni_agd) {
    vector[nint * max(int_id)] eta_agd = X_agd * beta_tilde;
    
    if (nint > 1) {  // If integration points are used
      // Loop over aux parameters first (i.e. by study) for efficiency
      for (i in 1:n_aux) {
        int ni = ni_aux_id_agd[i];
  
        if (ni) {
          array[ni] int wi = wi_aux_id_agd[i, 1:ni];
          
          for (j in 1:ni) {  // loop over times
            vector[nint] log_L_ii;

            log_L_ii = loglik_a(agd_time[wi[j]],
                                agd_itime[wi[j]],
                                agd_status[wi[j]],
                                eta_agd[((int_id[wi[j]]-1)*nint + 1):((int_id[wi[j]]-1)*nint + nint)],
                                scoef[i]);
  
            log_L_agd[wi[j]] = log_sum_exp(log_L_ii) - log(nint);
          }
        }
      }
    } else {  // If no integration
      // Loop over aux parameters (i.e. by study)
      for (i in 1:n_aux) {
        int ni = ni_aux_id_agd[i];
    
        if (ni) {
          array[ni] int wi = wi_aux_id_agd[i, 1:ni];
          log_L_agd[wi] = loglik(agd_time[wi],
                                 agd_itime[wi],
                                 agd_status[wi],
                                 eta_agd[int_id[wi]],
                                 scoef[i]);
        }
      }
    }
  }
}
model {
  // -- Priors --
  // These prior statements will cause Stan to raise warnings regarding transformed
  // parameters possibly needing Jacobian adjustments. These should be ignored, as
  // the transformation is entirely linear.
  mu ~ normal(0, prior_intercept_sd);
  beta ~ normal(0, prior_reg_sd);
  gamma ~ normal(0, prior_trt_sd);
  
  // RW1 prior on spline coefficients
  for (i in 1:n_aux) u_aux[i] ~ std_normal();
  sigma ~ normal(0, prior_sigma_sd);

  // -- Likelihood --
  target += log_L_ipd;
  target += log_L_agd;
}
generated quantities {
  // Log likelihood
  vector[ni_ipd + ni_agd] log_lik;

  if (ni_ipd) log_lik[1:ni_ipd] = log_L_ipd;
  if (ni_agd) log_lik[(ni_ipd + 1):(ni_ipd + ni_agd)] = log_L_agd;
}
