#' fit linear model to saliva, neutrophil and lymphocyte data for more than one participant
#' 
#' @param ids participant ids
#' @param neutrophils neutrophil population (same for all ids)
#' @param cells_lic.status (optional) list of character vectors specifying the cells and
#' licensing statuses to include for each id.  Each element should be a subset of
#' get_id_data(id)$cells_lic.status
#' @param model either a stan model object used to fit the model to data,
#' or a string specifying model name to compile this model
#' @param covariates character vector including "iKIR_count", "lic.status", "cells",
#' "disease"; or a single value "all" or "null"
#' @return list of two objects:
#' fit: stanfit object
#' data: data list used to fit stanfit object
#' @import dplyr
#' @export
fit_linear_stan <- function(ids, neutrophils, cells_lic.status, model, covariates,
                            adapt_delta = 0.8, seed = 1, init = "random") {
  
  stopifnot(all(covariates %in% c("cells", "lic.status", "iKIR_count", "disease",
                                  "all", "null")))
  # if using "all" covariates, list these out individually
  if(length(covariates) == 1 && covariates == "all") {
    covariates <- c("cells", "lic.status", "iKIR_count", "disease")
  }
  
  N <- length(ids)
  # if cell populations not specified, use all cell populations for each participant id
  if(missing(cells_lic.status)) {
    cells_lic.status <- lapply(ids, function(x) get_id_data(x)$cells_lic.status)
  }
  
  # read in saliva, granulocyte/monocyte and lympcyte data for each participant id
  extract_data_lists_for_id <- function(id, cells_lic.status) {
    data_in <- list(U = read_saliva_data(id),
                    Lb = read_gm_data(id, neutrophils),
                    L = read_lymphocyte_data(id))
    # exclude the t = 0 timepoint because it's 0 by definition of baseline correction -- doesn't help with fitting
    data_in <- lapply(data_in, function(x) filter(x, time > 0))
    
    # format the lymphocyte data
    pop <- length(cells_lic.status) # number of lymphocyte populations
    
    get_data_each_pop <- function(cells_lic.status1) {
      data_pop <- data_in$L %>%
        filter(cells_lic.status == cells_lic.status1)
      data_list <- list(  
        T_L= nrow(data_pop),
        t_L = data_pop$time,
        L = data_pop$L) # if 1, sample from posterior, if 0, sample from prior
    }
    
    data_pop <- lapply(cells_lic.status, get_data_each_pop)
    T_L <- vnapply(data_pop, function(x) x$T_L)
    T_Lmax <- max(T_L)
    T_Lsum <- sum(T_L)
    
    # format the saliva and granulocyte/monocyte data
    compartments <- c("U", "Lb")
    
    T_list <- lapply(compartments, function(x) nrow(data_in[[x]]))
    names(T_list) <- paste0("T_", compartments)
    t_list <- lapply(compartments, function(x) data_in[[x]]$time)
    names(t_list) <- paste0("t_", compartments)
    label_list <- lapply(compartments, function(x) data_in[[x]][[x]])
    names(label_list) <- compartments
    # define data and constants
    data_list <- c(T_list, t_list, label_list, 
                   list(data_pop = data_pop,
                        pop = pop,
                        T_Lmax = T_Lmax,
                        T_Lsum = T_Lsum)) 
  }
  
  # get the data for each participant id
  data_list <- Map(extract_data_lists_for_id, ids, cells_lic.status)
  
  # more formatting of data -- putting the data from each id into a single vector/matrix
  T_U <- vnapply(data_list, function(x) x$T_U)
  T_Umax <- max(T_U)
  T_Usum <- sum(T_U)
  T_Lb <- vnapply(data_list, function(x) x$T_Lb)
  T_Lbmax <- max(T_Lb)
  T_Lbsum <- sum(T_Lb)
  pop <- vnapply(data_list, function(x) x$pop)
  pop_max <- max(pop)
  C <- sum(pop)
  T_L <- numeric(C)
  U <- t_U <- matrix(0, nrow = T_Umax, ncol = N)
  Lb <- t_Lb <- matrix(0, nrow = T_Lbmax, ncol = N)
  j <- 1
  for(n in seq_len(N)) {
    t_U[seq_len(T_U[n]),n] <- data_list[[n]]$t_U
    U[seq_len(T_U[n]),n] <- data_list[[n]]$U
    t_Lb[seq_len(T_Lb[n]),n] <- data_list[[n]]$t_Lb
    Lb[seq_len(T_Lb[n]),n] <- data_list[[n]]$Lb
    for(p in seq_len(pop[n])) {
      T_L[j] <- data_list[[n]]$data_pop[[p]]$T_L
      j <- j + 1
    }
  }
  T_Lmax <- max(T_L)
  T_Lsum <- sum(T_L)
  t_L <- L <- array(0, dim = c(T_Lmax, C))
  C_to_N <- numeric(C)
  j <- 1
  for(n in seq_len(N)) {
    for(p in seq_len(pop[n])) {
      t_L[seq_len(T_L[j]),j] <- data_list[[n]]$data_pop[[p]]$t_L
      L[seq_len(T_L[j]),j] <- data_list[[n]]$data_pop[[p]]$L
      C_to_N[j] <- n
      j <- j + 1
    }
  }
  
  # define data and constants
  data_list <- list(  
    N = N,
    C = C,
    C_to_N = C_to_N,
    T_U = T_U,
    T_Lb = T_Lb,
    T_L = T_L,
    T_Umax = T_Umax,
    T_Lbmax = T_Lbmax,
    T_Lmax = T_Lmax,
    T_Usum = T_Usum,
    T_Lbsum = T_Lbsum,
    T_Lsum = T_Lsum,
    t_U = t_U,
    t_Lb = t_Lb,
    t_L = t_L,
    U = U,
    Lb = Lb,
    L = L)
  
  # read in values for fixed parameters for granulocytes/monocytes -- R and the delay from 
  # mitotic pool to blood
  fixed_pars <- get_fixed_pars(neutrophils)
  data_list <- c(data_list, as.list(fixed_pars)) 
  
  # read in functional iKIR count, licensing status, cell population and disease status
  # data if used as covariates
  if("iKIR_count" %in% covariates) {
    functional_iKIR_count <- lapply(ids, get_id_data) %>%
      vnapply(function(x) x[1,"functional_iKIR_count"])
    data_list <- c(data_list, list(iKIR_count = functional_iKIR_count))
  }
  if("lic.status" %in% covariates) {
    lic.status <- unlist(cells_lic.status)
    data_list <- c(data_list, list(is_unlicensed = as.numeric(grepl("Unlicensed", lic.status)),
                                   is_KIR_negative = as.numeric(grepl("KIR Negative", lic.status))))
  }
  if("cells" %in% covariates) {
    cells <- unlist(cells_lic.status)
    data_list <- c(data_list, list(is_TCM = as.numeric(grepl("TCM", cells))))
  }
  if("disease" %in% covariates) {
    disease <- lapply(ids, get_id_data) %>%
      vcapply(function(x) x[1,"disease"])
    data_list <- c(data_list, list(is_HCV = as.numeric(grepl("HCV", disease)),
                                   is_HIV = as.numeric(grepl("HIV", disease)),
                                   is_HTLV = as.numeric(grepl("HTLV", disease))))
  }
  # fit the model to data
  fit_stan_defaults(model, data_list, adapt_delta = adapt_delta, seed = seed, init = init)
}

#' predict all saliva, neutrophil and lymphocyte data for multiple participants
#' 
#' predict saliva and gm curves from model fit,
#' when saliva data and gm data for one cell population are fitted
#' together for a single participant
#' 
#' @param fit object returned by fit_saliva_gm_stan: list with two elements:
#' fit: stanfit object
#' data: data list used to fit stanfit object
#' @param saliva_model either a stan model object used to predict saliva curves,
#' or a string specifying model name to compile this model
#' @param gm_model either a stan model object used to predict gm curves,
#' or a string specifying model name to compile this model
#' @param N_samples integer.  number of samples to use in model predictions
#' @param lymphocyte_model either a stan model object used to predict lymphocyte curves,
#' or a string specifying model name to compile this model
#' @return list of three objects:
#' time: time series for which to predict
#' N_predictions: list of length N, where each element is for one participants.
#' Each element is in turn a list of two elements:
#' Uhat: matrix of model predictions.  Each row is a time and each column is a sample.
#' Lbhat: matrix of model predictions.  Each row is a time and each column is a sample.
#' C_predictions: list of length C, where each element is for one cell population.
#' Each element is a matrix of model predictions.  Each row is a time and each column is a sample.  
#' @importFrom dplyr %>%
#' @export
pred_null_from_fit <- function(fit, saliva_model, gm_model, lymphocyte_model, N_samples = 100) {
  # read in N_samples samples from the posterior distribution (if N_samples is greater than the number of samples
  # in the posterior, take the minimum of the two)
  ext_fit <- extract_fit(fit$fit, N_samples = N_samples)
  # adjust N_samples to reflect the actual number of samples taken
  N_samples <- nrow(ext_fit)
  # predict every day until the latest data point across all individuals and compartments
  pred_times <- seq(0, max(c(100, fit$data$t_U, fit$data$t_L, fit$data$t_Lb)))
  
  # get the parameter names associated with each compartment (saliva, granulocytes/monocytes, lymphocytes)
  par_names_compartment <- get_par_names()
  # number of participants
  N <- fit$data$N
  # number of lymphocyte populations
  C <- fit$data$C
  
  # fitted parameter names which are individual-specific
  scalar_par_names <- "delta"
  N_par_names <- c(par_names_compartment$U, 
                   par_names_compartment$Lb)
  N_par_names <- N_par_names[!(N_par_names %in% scalar_par_names)]
  
  # fitted parameter names which are lymphocyte population specific
  C_par_names <- par_names_compartment$L
  if("p_sd" %in% colnames(ext_fit)) {
    C_par_names[C_par_names == "pb_w"] <- "p"
  }
  
  # fixed parameter names for granulocytes/monocytes
  constant_names <- get_fixed_pars("Granulocytes") %>% names
  
  # gather parameter values for granulocyte/monocyte label predictions
  make_data_list_gm <- function(n) {
    data_list_scalar_pars <- lapply(scalar_par_names, function(x) array(ext_fit[, x]))
    data_list_vector_pars <- lapply(N_par_names, 
                                    function(x) array(ext_fit[, index_par_name(x, n)]))
    names(data_list_scalar_pars) <- scalar_par_names
    names(data_list_vector_pars) <- N_par_names
    data_list_constants <- fit$data[constant_names]
    names(data_list_constants) <- constant_names
    data_list_gm <- c(data_list_scalar_pars, data_list_vector_pars, data_list_constants, 
                      list(N_samples = N_samples,
                           T= length(pred_times),
                           ts = pred_times))
    data_list_gm
  }
  
  # gather parameter values for saliva and granulocyte/monocyte label predictions
  make_data_list_n <- function(n) {
    data_list_gm <- make_data_list_gm(n)
    gm_only_pars <- c("z", "b_w", "R", "deltaT")
    
    data_list_saliva <- data_list_gm[!(names(data_list_gm) %in% gm_only_pars)]
    list(U = data_list_saliva, Lb = data_list_gm)
  }
  
  # solve for fraction of label in saliva and granulocytes/monocytes
  make_N_predictions <- function(n) {
    
    data_list <- make_data_list_n(n)
    
    Uhat <- make_pred(saliva_model, data_list$U) %>%
      extract_pred(N_samples = N_samples)
    Lbhat <- make_pred(gm_model, data_list$Lb) %>%
      extract_pred(N_samples = N_samples)
    
    list(Uhat = Uhat, Lbhat = Lbhat)
  }
  
  # solve for fraction of label in lymphocytes
  make_C_predictions <- function(c1) {
    n <- fit$data$C_to_N[c1]
    data_list_n <- make_data_list_n(n)
    data_list_saliva <- data_list_n$U
    data_list_lymphocyte_pop <- lapply(C_par_names,
                                       function(x) array(ext_fit[, index_par_name(x, c1)]))
    names(data_list_lymphocyte_pop) <- C_par_names
    
    
    data_list_lymphocyte_pop$pb_w <- data_list_lymphocyte_pop$p * data_list_n$Lb$b_w
    
    data_list_lymphocyte_pop <- data_list_lymphocyte_pop[names(data_list_lymphocyte_pop) != "p"]
    data_list <- c(data_list_saliva, data_list_lymphocyte_pop)
    
    Lhat <- make_pred(lymphocyte_model, data_list) %>%
      extract_pred(N_samples = N_samples)
    
    Lhat
  }
  
  N_predictions <- lapply(seq_len(N), make_N_predictions)
  C_predictions <- lapply(seq_len(C), make_C_predictions)
  
  pred <- list(time = pred_times,
               N_predictions = N_predictions,
               C_predictions = C_predictions)
  pred
}