// function to calculate the fraction of label in saliva over time
functions {
  // declaration necessary before defining contents of function
  // because it's a recursive function
      real calc_U(real t,
             real frac,
             real delta,
             real phase_one_end,
             real label_end);
             
    real calc_U(real t,
               real frac,
               real delta,
               real phase_one_end,
               real label_end) {
      real U;
      real t0;
      real f1;
      // returns analytic solution to fraction of label in body water
      if(t == 0) {
        U = 0;
        return U;
      } else if(t <= phase_one_end) {
        f1 = frac;
        t0 = 0;
      } else if (t <= label_end) {
        // in the protocol, the amount of heavy water is reduced by 1/3 after the first week
        f1 = 2./3. * frac;
        t0 = phase_one_end;
      } else {
        f1 = 0;
        t0 = label_end;
      }
      U = f1 * (1 - exp(-delta * (t - t0))) + calc_U(t0, frac, delta, phase_one_end, label_end) * exp(-delta * (t - t0));
      return U;
    }
    
    // calculate the fraction of label in granulocytes/monocytes in the mitotic pool (i.e. without delay
    // to being observed in the blood)
    real calc_X(real t,
     real frac,
     real delta,
     real phase_one_end,
     real label_end,
     real z,
     real R,
     real b_w);
    real calc_X(real t,
             real frac,
             real delta,
             real phase_one_end,
             real label_end,
             real z,
             real R,
             real b_w) {
    real X;
    real A;
    real B;
    real C;
    real D;
    real E;
    real F;
    real G;
    real H;
    real J;
    real L_p_phase_one_end;
    real L_p_label_end;
    real t_shift;
    real t_inter;
    A = b_w * frac;
    B = A * delta / (z * R - delta);
    J = 2. /3. * frac;
    H = calc_U(phase_one_end, frac, delta, phase_one_end, label_end) - J;
    L_p_phase_one_end = A + B * exp(-z * R * phase_one_end) - (A + B) * exp(-delta * phase_one_end);
    C = z * R * b_w * H / (z * R - delta);
    E = J * b_w;
    D = -(C + E) + L_p_phase_one_end;
    if(t <= phase_one_end) {
      X = A * (1 - exp(-z * t)) + B / (1 - R) * (exp(-z * R * t) - exp(-z * t)) - 
      (A + B) * z / (z - delta) * (exp(-delta * t) - exp(-z * t));
    } else if (t <= label_end) {
      t_shift = t - phase_one_end;
      X = calc_X(phase_one_end, frac, delta, phase_one_end, label_end, z, R, b_w) * exp(-z * t_shift) +
      z * C / (z - delta) * (exp(-delta * t_shift) - exp(-z * t_shift)) +
      D / (1 - R) * (exp(-z * R * t_shift) - exp(-z * t_shift)) +
      E * (1 - exp(-z * t_shift));
    } else {
      t_inter = label_end - phase_one_end;
      L_p_label_end = C * exp(-delta * t_inter) + D * exp(-z * R * t_inter) + E;
      G = z * R * b_w * calc_U(label_end, frac, delta, phase_one_end, label_end) / (z * R - delta);
      F = L_p_label_end - G;
      t_shift = t - label_end;
      X = calc_X(label_end, frac, delta, phase_one_end, label_end, z, R, b_w) * exp(-z * t_shift) + 
      F / (1 - R) * (exp(-z * R * t_shift) - exp(-z * t_shift)) +
      z * G / (z - delta) * (exp(-delta * t_shift) - exp(-z * t_shift));
    }
    return X;
  }
  
  // calculate fraction of label in granulocytes/monocytes in the blood -- i.e. add the delay onto
  // the fraction of label in the mitotic pool
      real calc_Lb(real t,
           real frac,
           real delta,
           real phase_one_end,
           real label_end,
           real z,
           real R,
           real b_w,
           real deltaT);
           
    real calc_Lb(real t,
           real frac,
           real delta,
           real phase_one_end,
           real label_end,
           real z,
           real R,
           real b_w,
           real deltaT) {
             real Lb;
             if(t <= deltaT) {
               Lb = 0;
             } else {
                Lb = calc_X(t - deltaT, frac, delta, phase_one_end, label_end, z, R, b_w);
             }

           return Lb;
           }
  
  // calculate fraction of label in lymphocytes in the lymph nodes
    real calc_Y(real t,
       real frac,
       real delta,
       real phase_one_end,
       real label_end,
       real pb_w,
       real dstar);
    real calc_Y(real t,
             real frac,
             real delta,
             real phase_one_end,
             real label_end,
             real pb_w,
             real dstar) {
    real Y;
    real A;
    real B;
    real t_shift;
    if(t <= phase_one_end) {
      Y = pb_w * frac / (delta - dstar) * (delta/dstar *(1 - exp(-dstar * t)) - (1 - exp(-delta * t)));
    } else if (t <= label_end) {
      B = 2. /3. * frac;
      A = calc_U(phase_one_end, frac, delta, phase_one_end, label_end) - B;
      t_shift = t - phase_one_end;
      Y = calc_Y(phase_one_end, frac, delta, phase_one_end, label_end, pb_w, dstar) * exp(-dstar * t_shift) +
      pb_w * exp(-dstar * t_shift) * (A / (dstar - delta) *(exp((dstar - delta) * t_shift) - 1) +
      B / dstar * (exp(dstar * t_shift) - 1));
    } else {
      t_shift = t - label_end;
      Y = calc_Y(label_end, frac, delta, phase_one_end, label_end, pb_w, dstar) * exp(-dstar * t_shift) + pb_w * exp(-dstar * t_shift) * calc_U(label_end, frac, delta, phase_one_end, label_end) / (dstar - delta) * (exp((dstar - delta) * t_shift) - 1);
    }
    return Y;
  }
  
  // calculate fraction of label in lymphocytes in the blood
  real calc_L(real t,
         real frac,
         real delta,
         real phase_one_end,
         real label_end,
         real pb_w,
         real dstar,
         real delay) {
           real L;
           if(t <= delay) {
             L = 0;
           } else {
              L = calc_Y(t - delay, frac, delta, phase_one_end, label_end, pb_w, dstar);
           }

         return L;
         }

}

// data and fixed parameter values
data {
  int<lower=1> N; // number of participants
  int<lower=N> C; // total number of cell populations across participants
  int C_to_N[C]; // which N each cell pop belongs to
  int T_U[N]; // number of sampling times for saliva
  int T_Lb[N]; // number of sampling times for neutrophils
  int T_L[C]; // number of sampling times for lymphocytes
  int<lower=1> T_Umax; // maximum number of sampling times for saliva across participants
  int<lower=1> T_Lbmax; // maximum number of sampling times for neutrophils across participants
  int<lower=1> T_Lmax; // maximum number of sampling times for lymphocytes across participants
  int<lower=1> T_Usum; // total number of sampling times across saliva
  int<lower=1> T_Lbsum; // total number of sampling times across neutrophils
  int<lower=1> T_Lsum; // total number of sampling times across lymphocytes
  matrix[T_Umax, N] t_U; // sampling times
  matrix[T_Lbmax, N] t_Lb; // sampling times
  matrix[T_Lmax, C] t_L; // sampling times
  matrix[T_Umax, N] U; // fraction of label in saliva
  matrix[T_Lbmax, N] Lb; // fraction of label in neutrophils
  matrix[T_Lmax, C] L; // fraction of label in granulocytes
  real<lower=0> phase_one_end; // time at which label intake decreases in the protocol
  real<lower=phase_one_end> label_end; // time at which labelling ends
  real<lower = 0.> R;// ratio of mitotic pool to blood granulocyte/monocyte population size
  real<lower = 0.> deltaT; // delay between fraction of label in mitotic pool and that in blood for
  // granulocytes/monocytes
  int evaluate_likelihood; // if 1, sample from posterior, if 0, sample from prior
}

// fitted parameters
parameters {
  // ranges are ranges of uniform priors
  vector<lower=0.,upper=0.1>[N] frac; // maximum possible fraction of label in body water
  real<lower=0.,upper=1.> delta; // turnover rate of body water
  vector<lower = 0., upper = 10.>[N] z; // turnover rate of granulocytes/monocytes
  vector<lower = 0., upper = 7.>[N] b_w;
  real<upper = log(.05)> p_baseline; // mu parameter in lognormal distribution for p (hyperparameter)
  real<upper = log(1.)> dstar_mean; // mu parameter in lognormal distribution for dstar  (hyperparameter)
  real<upper = log(21.)> delay_mean; // mu parameter in lognormal distribution for delay for lymphocytes (hyperparameter)
  real<lower = 0., upper = 1.> p_sd; // sigma parameter in lognormal distribution for p (hyperparameter)
  real<lower = 0., upper = 1.> dstar_sd; // sigma parameter in lognormal distribution for dstar (hyperparameter)
  real<lower = 0., upper = 1.> delay_sd; // sigma parameter in lognormal distribution for delay (hyperparameter)
  real<lower=-5.,upper = -1.> log10_sigma_U; // log10 of the standard deviation of observation model for saliva
  real<lower=-5.,upper = -1.> log10_sigma_Lb; // log10 of the standard deviation of observation model for granulocytes/monocytes
  real<lower=-5.,upper = -1.> log10_sigma_L; // log10 of the standard deviation of observation model for lymphocytes
  vector<lower=0.>[C] p; // p for each of the cell populations
  vector<lower=0.>[C] dstar; // dstar for each of the cell populations
  vector<lower=0.>[C] delay; // delay for each of the cell populations
}

transformed parameters {
  vector[C] pb_w; // only p * b_w appears in the model equations -- precalculate this
  matrix[T_Umax, N] Uhat; // solution to fraction of label in saliva for sampled parameters
  matrix[T_Lbmax, N] Lbhat; // solution to fraction of label in granulocytes/monocytes for sampled parameters
  matrix[T_Lmax, C] Lhat; // solution to fraction of label in lymphocytes for sampled parameters

  real sigma_U; // standard deviation  of observation model for saliva
  real sigma_Lb; // standard deviation  of observation model for granulocytes/monocytes
  real sigma_L; // standard deviation  of observation model for lymphocytes
  
  sigma_U = 10.^log10_sigma_U;
  sigma_Lb = 10.^log10_sigma_Lb;
  sigma_L = 10.^log10_sigma_L;

  pb_w = p .* b_w[C_to_N];

  Uhat = rep_matrix(0., T_Umax, N);
  Lbhat = rep_matrix(0., T_Lbmax, N);
  Lhat = rep_matrix(0., T_Lmax, C);

  // calculate solution to fraction of label for sampled parameters
  if(evaluate_likelihood) {
    for(n in 1:N) {
      for (t in 1:T_U[n]) {
        Uhat[t,n] = calc_U(t_U[t,n], frac[n], delta, phase_one_end, label_end);
      }
      for (t in 1:T_Lb[n]) {
        Lbhat[t,n] = calc_Lb(t_Lb[t,n], frac[n], delta, phase_one_end, label_end, z[n], R, b_w[n], deltaT);
      }
    }
    for(c in 1:C) {
      for (t in 1:T_L[c]) {
        Lhat[t,c] = calc_L(t_L[t,c], frac[C_to_N[c]], delta, phase_one_end, label_end, pb_w[c], dstar[c], delay[c]);
      }
    }

  }
}

// calculate likelihood for model parameters
model {
  // p, dstar and delay are lognormally distributed; the parameters of teh lognormal distributions are fitted
  // hypoerparameters
  p ~ lognormal(p_baseline, p_sd);
  dstar ~ lognormal(dstar_mean, dstar_sd);
  delay ~ lognormal(delay_mean, delay_sd);

  // the fraction of label is normally distributed around the model solution
  if(evaluate_likelihood) {
    // loop over individuals for saliva and granulocytes/monocytes
    for(n in 1:N) {
      for (t in 1:T_U[n]) {
        U[t,n] ~ normal(Uhat[t,n], sigma_U);
      }
      for (t in 1:T_Lb[n]) {
        Lb[t,n] ~ normal(Lbhat[t,n], sigma_Lb);
      }
    }
    // loop over cell populations/licensing statuses/individuals for lymphocytes
      for(c in 1:C) {
        for (t in 1:T_L[c]) {
          L[t,c] ~ normal(Lhat[t,c], sigma_L);
        }
      }
  }
}

// because of the way Stan works, if we want to save the likelihood at each step
// we need to recalculate it
generated quantities {
  vector[T_Usum + T_Lbsum + T_Lsum] log_lik;
  int j;
  j = 1;

  for(n in 1:N) {
    for (t in 1:T_U[n]) {
      log_lik[j] = normal_lpdf(U[t,n] | Uhat[t,n], sigma_U);
      j = j + 1;
    }
    for (t in 1:T_Lb[n]) {
      log_lik[j] = normal_lpdf(Lb[t,n] | Lbhat[t,n], sigma_Lb);
      j = j + 1;
    }
  }
  
    for(c in 1:C) {
      for (t in 1:T_L[c]) {
        log_lik[j] = normal_lpdf(L[t,c] | Lhat[t,c], sigma_L);
        j = j + 1;
      }
    }
}
