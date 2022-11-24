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
        // in the protocol, the amount of heavy water is reduced by 1/3 after the first week
      } else if (t <= label_end) {
        f1 = 2./3. * frac;
        t0 = phase_one_end;
      } else {
        f1 = 0;
        t0 = label_end;
      }
      U = f1 * (1 - exp(-delta * (t - t0))) + calc_U(t0, frac, delta, phase_one_end, label_end) * exp(-delta * (t - t0));
      return U;
    }
  // calculate fraction of label in lymphocytes in the lymph nodes
    real calc_X(real t,
         real frac,
         real delta,
         real phase_one_end,
         real label_end,
         real pb_w,
         real dstar);
    real calc_X(real t,
             real frac,
             real delta,
             real phase_one_end,
             real label_end,
             real pb_w,
             real dstar) {
    real X;
    real A;
    real B;
    real t_shift;
    if(t <= phase_one_end) {
      X = pb_w * frac / (delta - dstar) * (delta/dstar *(1 - exp(-dstar * t)) - (1 - exp(-delta * t)));

    } else if (t <= label_end) {
      B = 2. /3. * frac;
      A = calc_U(phase_one_end, frac, delta, phase_one_end, label_end) - B;
      t_shift = t - phase_one_end;
      X = calc_X(phase_one_end, frac, delta, phase_one_end, label_end, pb_w, dstar) * exp(-dstar * t_shift) +
      pb_w * exp(-dstar * t_shift) * (A / (dstar - delta) *(exp((dstar - delta) * t_shift) - 1) +
      B / dstar * (exp(dstar * t_shift) - 1));
    } else {
      t_shift = t - label_end;
      X = calc_X(label_end, frac, delta, phase_one_end, label_end, pb_w, dstar) * exp(-dstar * t_shift) + pb_w * exp(-dstar * t_shift) * calc_U(label_end, frac, delta, phase_one_end, label_end) / (dstar - delta) * (exp((dstar - delta) * t_shift) - 1);
    }
    return X;
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
                L = calc_X(t - delay, frac, delta, phase_one_end, label_end, pb_w, dstar);
             }

           return L;
           }
}

// data and fixed parameter values
data {
  int<lower=1> N_samples; // number of posterior samples
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real<lower=0> phase_one_end;// time at which label intake decreases in the protocol
  real<lower=phase_one_end> label_end; // time at which labelling ends
  real<lower=0.> frac[N_samples]; // maximum possible fraction of label in body water
  real<lower=0.> delta[N_samples]; // turnover rate of body water
  real<lower = 0.> pb_w[N_samples]; // average proliferation rate of lymphocytes * b_w
  real<lower = 0.> dstar[N_samples]; // disappearance rate of fast compartment
  real<lower = 0.> delay[N_samples]; // delay between uptake of label in lymph node and appearance of label in blood
}

// model parameters block is empty because we're not performing inference
parameters {
}

// model block is empty because we're not performing inference
model {
}

// calculate fraction of label
generated quantities {
  real Lhat[T,N_samples];
  for (t in 1:T) {
    for (n in 1:N_samples) {
        Lhat[t,n] = calc_L(ts[t], frac[n], delta[n], phase_one_end, label_end, pb_w[n], dstar[n], delay[n]);
    }
  }
}
