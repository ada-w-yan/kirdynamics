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
}

// data and fixed parameter values
data {
  int<lower=1> N_samples; // number of samples from posterior for which to calculate fraction of label
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real<lower=0> phase_one_end; // time at which label intake decreases in the protocol
  real<lower=phase_one_end> label_end; // time at which labelling ends
  real<lower=0.,upper=0.1> frac[N_samples]; // maximum possible fraction of label in body water
  real<lower=0.,upper=1.> delta[N_samples]; // turnover rate of body water
}

// model parameters block is empty because we're not performing inference
parameters {
}
// model block is empty because we're not performing inference
model {
}

// calculate fraction of label
generated quantities {
  real Uhat[T,N_samples];
  for (t in 1:T) {
    for (n in 1:N_samples) {
        Uhat[t,n] = calc_U(ts[t], frac[n], delta[n], phase_one_end, label_end);
    }
  }
}
