// function to calculate the fraction of label in saliva over time
functions {
    // declaration necessary before defining contents of function
  // because it's a recursive function
      real calc_U(real t,
             real frac,
             real delta,
             real phase_one_end,
             real label_end);
      // returns analytic solution to fraction of label in body water           
    real calc_U(real t,
               real frac,
               real delta,
               real phase_one_end,
               real label_end) {
      real U;
      real t0;
      real f1;
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

}

// data and fixed parameter values
data {
  int<lower=1> N_samples; // number of samples from posterior for which to calculate fraction of label
  int<lower=1> T; // number of sampling times
  real ts[T]; // sampling times
  real<lower=0> phase_one_end; // time at which label intake decreases in the protocol
  real<lower=phase_one_end> label_end; // time at which labelling ends
  real<lower=0.,upper=0.1> frac[N_samples]; // maximum possible fraction of label in body water
  real<lower=0.> delta[N_samples]; // turnover rate of body water
  real<lower = 0.> R; // ratio of mitotic pool to blood granulocyte/monocyte population size
  real<lower = 0.> deltaT;  // delay between fraction of label in mitotic pool and that in blood for
  // granulocytes/monocytes
  real<lower = 0.> z[N_samples]; // turnover rate of granulocytes/monocytes
  real<lower = 0., upper = 7.> b_w[N_samples];
}

// model parameters block is empty because we're not performing inference
parameters {
}
// model block is empty because we're not performing inference
model {
}

// calculate fraction of label
generated quantities {
  real Lbhat[T,N_samples];
  for (t in 1:T) {
    for (n in 1:N_samples) {
        Lbhat[t,n] = calc_Lb(ts[t], frac[n], delta[n], phase_one_end, label_end, z[n], R, b_w[n], deltaT);
    }
  }
}
