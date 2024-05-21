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
  
  // -- Generalised Gamma ldpf --
  real gengamma_lpdf(real y, real mu, real sigma, real k) {
    // Parameterisation of Lawless
	  real Q = pow(k, -0.5);
	  real w = Q * (log(y) - mu) / sigma;
	  return -log(sigma) - log(y) - 0.5 * log(k) * (1 - 2 * k) + k * (w - exp(w)) - lgamma(k);
  }

  //-- log Survival functions --
  // For efficiency we need two forms: one vectorised over individuals (IPD), 
  // one vectorised over integration points (AgD)
  
  vector lS_ipd(int dist, vector y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    int n = num_elements(y);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = -y .* exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)./aux .* expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -y .* exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) .* exp(-aux * eta);
    } else if (dist == 6) { // log Normal
      for (i in 1:n) out[i] = lognormal_lccdf(y[i] | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y[i] | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      real Q = inv(sqrt(aux2));
      vector[n] w = exp(Q * (log(y) - eta) / aux) * aux2;
      out = log1m(gamma_p(aux2, w));
    }

    return out;
  }
  
  vector lS_agd(int dist, real y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k
    
    int n = num_elements(eta);
    vector[n] out;
    
    if (dist == 1) { // Exponential
      out = -y * exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow(y, aux) * exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)/aux * expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -y * exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) * exp(-aux * eta);
    } else if (dist == 6) { // log Normal
      for (i in 1:n) out[i] = lognormal_lccdf(y | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      real Q = inv(sqrt(aux2));
      vector[n] w = exp(Q * (log(y) - eta) / aux) * aux2;
      out = log1m(gamma_p(aux2, w));
    }

    return out;
  }
  
  // -- log Hazard functions --
  // For efficiency we need two forms: one vectorised over individuals (IPD), 
  // one vectorised over integration points (AgD)
  
  vector lh_ipd(int dist, vector y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    int n = num_elements(y);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = eta;
    } else if (dist == 2) { // Weibull
      out = log(aux) + eta + lmultiply(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      out = eta + (aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -eta;
    } else if (dist == 5) { // Weibull AFT
      out = log(aux) - (aux * eta) + lmultiply(aux - 1, y);
    } else if (dist == 6) { // log Normal
      for (i in 1:n) out[i] = lognormal_lpdf(y[i] | eta[i], aux) - lognormal_lccdf(y[i] | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y[i] | aux, eeta[i]) - gamma_lccdf(y[i] | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }
  
  vector lh_agd(int dist, real y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape
    
    int n = num_elements(eta);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = eta;
    } else if (dist == 2) { // Weibull
      out = log(aux) + eta + lmultiply(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      out = eta + (aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -eta;
    } else if (dist == 5) { // Weibull AFT
      out = log(aux) - (aux * eta) + lmultiply(aux - 1, y);
    } else if (dist == 6) { // log Normal
      for (i in 1:n) out[i] = lognormal_lpdf(y | eta[i], aux) - lognormal_lccdf(y | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y | aux, eeta[i]) - gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }
  
  // -- log likelihood function --
  // For efficiency we need two forms: one vectorised over individuals (IPD), 
  // one vectorised over integration points (AgD)
  
  vector loglik_ipd(int dist, vector time, array[] int status, vector eta, real aux, real aux2) {
    int n = num_elements(eta);
    vector[n] l;

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      array[n - sum(status)] int w0 = which(status, 0);
      
      // Censored
      l[w0] = lS_ipd(dist, time[w0], eta[w0], aux, aux2);

      // Observed events
      if (dist == 6) for (i in 1:n) if (status[i] == 1) l[i] = lognormal_lpdf(time[i] | eta[i], aux);
      if (dist == 8) {
        vector[n] eeta = exp(-eta);
        for (i in 1:n) if (status[i] == 1) l[i] = gamma_lpdf(time[i] | aux, eeta[i]);
      }
      if (dist == 9) for (i in 1:n) if (status[i] == 1) l[i] = gengamma_lpdf(time[i] | eta[i], aux, aux2);

    } else {
      array[sum(status)] int w1 = which(status, 1);

      // All observations get survival contribution
      l = lS_ipd(dist, time, eta, aux, aux2);

      // Observed events get additional hazard contribution
      l[w1] += lh_ipd(dist, time[w1], eta[w1], aux, aux2);
    }

    return l;
  }
  
  vector loglik_agd(int dist, real time, int status, vector eta, real aux, real aux2) {
    
    int n = num_elements(eta);
    vector[n] l;

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      if (status == 0) { // Censored
        l = lS_agd(dist, time, eta, aux, aux2);
      } else { // Observed
        if (dist == 6) for (i in 1:n) l[i] = lognormal_lpdf(time | eta[i], aux);
        if (dist == 8) {
          vector[n] eeta = exp(-eta);
          for (i in 1:n) l[i] = gamma_lpdf(time | aux, eeta[i]);
        }
        if (dist == 9) for (i in 1:n) l[i] = gengamma_lpdf(time | eta[i], aux, aux2);
      }
    } else {
      if (status == 0) { // Censored
        l = lS_agd(dist, time, eta, aux, aux2);
      } else { // Observed
        l = lS_agd(dist, time, eta, aux, aux2) + lh_agd(dist, time, eta, aux, aux2);
      }
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
  int<lower=0> nX; // number of columns of design matrix
  
  // -- Survival data --
  vector<lower=0>[ni_ipd + ni_agd] time; // observation time
  array[ni_ipd + ni_agd] int<lower=0, upper=1> status; // event status (0 = censored, 1 = event)
  
  // -- Auxiliary parameter details --
  array[ni_ipd + ni_agd] int<lower=1> aux_id; // aux id (e.g. study id for a PH/AFT model)

  // -- Integration point details --
  array[ni_agd] int<lower=1> int_id; // integration sample id (i.e. for an arm or study)
  
  // -- Distribution flag --
  int<lower=1, upper=9> dist;
  // 1 = Exponential PH
  // 2 = Weibull PH
  // 3 = Gompertz PH
  // 4 = Exponential AFT
  // 5 = Weibull AFT
  // 6 = log Normal AFT
  // 7 = log logistic AFT
  // 8 = Gamma AFT
  // 9 = Generalised Gamma AFT

  // -- Design matrix --
  // Either X, or QR Q if QR=TRUE
  matrix[ni_ipd + (ni_agd ? nint * max(int_id) : 0), nX] X;
  
  // -- QR decomposition --
  int<lower=0, upper=1> QR; // QR decomposition flag (1=used QR decomposition)
  matrix[QR ? nX : 0, QR ? nX : 0] R_inv;
  
  // -- Priors --
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_reg_sd;
  real<lower=0> prior_trt_sd;
  real<lower=0> prior_aux_sd;
  real<lower=0> prior_aux2_sd;
}
transformed data {
  // Exponential model indicator, 0 = exponential
  int<lower=0, upper=1> nonexp = (dist == 1 || dist == 4) ? 0 : 1;
  // Generalised Gamma model indicator, 1 = gengamma
  int<lower=0, upper=1> gengamma = dist >= 9 ? 1 : 0;
  int n_aux_id = max(aux_id); // Number of auxiliary parameters (e.g. number of studies)

  // Split out IPD and AgD 
  vector[ni_ipd] time_ipd = head(time, ni_ipd);
  vector[ni_agd] time_agd = tail(time, ni_agd);
  
  array[ni_ipd] int status_ipd = head(status, ni_ipd);
  array[ni_agd] int status_agd = tail(status, ni_agd);
  
  matrix[ni_ipd, nX] X_ipd;
  matrix[ni_agd ? nint * max(int_id) : 0, nX] X_agd;
  
  array[ni_ipd] int aux_id_ipd = nonexp ? head(aux_id, ni_ipd) : rep_array(1, ni_ipd);
  array[ni_agd] int aux_id_agd = nonexp ? tail(aux_id, ni_agd) : rep_array(1, ni_agd);
  
  // -- Aux ID indexing arrays --
  // For efficiency we will loop over aux_id's, vectorising within each group

  // Number of observations within each ID
  array[n_aux_id] int ni_aux_id_ipd = nwhich_all(aux_id_ipd, n_aux_id);
  array[n_aux_id] int ni_aux_id_agd = nwhich_all(aux_id_agd, n_aux_id);
  
  // Index lookup - which observations belong to group i?
  array[n_aux_id, max(ni_aux_id_ipd)] int wi_aux_id_ipd;
  array[n_aux_id, max(ni_aux_id_agd)] int wi_aux_id_agd;

  for (i in 1:n_aux_id) {
    if (ni_aux_id_ipd[i]) wi_aux_id_ipd[i, 1:ni_aux_id_ipd[i]] = which(aux_id_ipd, i);
    if (ni_aux_id_agd[i]) wi_aux_id_agd[i, 1:ni_aux_id_agd[i]] = which(aux_id_agd, i);
  }
  
  // Split out IPD and AgD
  if (ni_ipd) X_ipd = X[1:ni_ipd, ];
  if (ni_agd) X_agd = X[(ni_ipd + 1):(ni_ipd + nint * max(int_id)), ];
}
parameters {
  // Parameters on QR scale
  vector[nX] beta_tilde;
  
  // Auxiliary parameters for parametric model (typically shape)
  // Exponential model has shape = 1 so parameter is removed (zero dim)
  vector<lower=0>[nonexp ? ns_ipd + ns_agd : 0] aux;
  // Second auxiliary parameter for generalised gamma
  vector<lower=0>[gengamma ? ns_ipd + ns_agd : 0] aux2;
}
transformed parameters {
  // -- Back-transformed parameters --
  vector[nX] allbeta = QR ? R_inv * beta_tilde : beta_tilde;
  // Study intercepts
  vector[ns_ipd + ns_agd] mu = allbeta[1:(ns_ipd + ns_agd)];
  // Treatment effects
  vector[nt - 1] gamma = allbeta[(ns_ipd + ns_agd + 1):(ns_ipd + ns_agd + nt - 1)];
  // Regression terms
  vector[nX - (nt-1) - ns_ipd - ns_agd] beta = allbeta[(ns_ipd + ns_agd + nt):nX];
  
  // -- log likelihood contributions --
  vector[ni_ipd + ni_agd] log_lik;
  
  // -- Evaluate IPD likelihood --
  if (ni_ipd) {
    vector[ni_ipd] eta_ipd = X_ipd * beta_tilde;  // Linear predictor
    
    // Loop over aux IDs
    for (i in 1:n_aux_id) {
        int ni = ni_aux_id_ipd[i];

        if (ni) {
          array[ni] int wi = wi_aux_id_ipd[i, 1:ni];
          real auxi;
          real aux2i;
          if (nonexp) auxi = aux[i];
          if (gengamma) aux2i = aux2[i];

          log_lik[wi] = loglik_ipd(dist,
                                   time_ipd[wi],
                                   status_ipd[wi],
                                   eta_ipd[wi],
                                   auxi,
                                   aux2i);
        }
      }
  }
  
  // -- Evaluate AgD likelihood --
  if (ni_agd) {
    vector[nint * max(int_id)] eta_agd = X_agd * beta_tilde;  // Linear predictor
    
    // Loop over aux IDs
    for (i in 1:n_aux_id) {
      int ni = ni_aux_id_agd[i];

      if (ni) {
        array[ni] int wi = wi_aux_id_agd[i, 1:ni];
        real auxi = nonexp ? aux[i] : 0;
        real aux2i = gengamma ? aux2[i] : 0;

        // Loop over observations, evaluating marginal likelihood for each via integration
        for (j in 1:ni) {
          vector[nint] log_L_ii;

          log_L_ii = loglik_agd(dist,
                                time_agd[wi[j]],
                                status_agd[wi[j]],
                                eta_agd[((int_id[wi[j]]-1)*nint + 1):((int_id[wi[j]]-1)*nint + nint)],
                                auxi,
                                aux2i);

          log_lik[ni_ipd + wi[j]] = log_sum_exp(log_L_ii) - log(nint);
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
  aux ~ normal(0, prior_aux_sd);
  aux2 ~ normal(0, prior_aux2_sd);
  
  // -- Likelihood --
  target += log_lik;
}
