library(ggplot2)
library(mrgsolve)
library(survival) # for tcc plot
library(survminer) # for tcc plot

# Lesion PK model (IV Dosing)
code <-"$PARAM
        ka = 1, 
        CL  = 1,
        V  = 10,
        PC  = 1,
        KPL = 10
      
        $CMT
          GUT
          CP
          LESION
        
        $ODE
          dxdt_GUT  = -ka * GUT;
          dxdt_CP   = ka * GUT -(CL/V)*CP;
          dxdt_LESION =  KPL * (PC*CP/V - LESION);
        
        $CAPTURE
          C_CP      = CP/V;
          C_LESION  = LESION;"

pkmod <- mcode("plasma_pk", code)

# Run single regimen PK simulation
run_pk_sim <- function(target_cmt, ka, clin_CL, clin_V, dose, ndoses, inter) {
  
  regimen <- ev(
    amt = dose, 
    ii = inter, 
    addl = ndoses - 1,
    cmt  = target_cmt)
  
  out <- pkmod %>%
    param(
      ka = ka, 
      CL = clin_CL, 
      V = clin_V) %>%
    ev(regimen) %>%
    mrgsim(end = inter * ndoses, delta = 0.1)
  
  as_tibble(out)
}



# Run single regimen lesion PK simulation
run_lesion_pk_sim <- function(target_cmt, ka, clin_CL, clin_V, dose, ndoses, inter, PC, KPL = 10) {
  
  regimen <- ev(
    amt = dose, 
    ii = inter, 
    addl = ndoses - 1,
    cmt  = target_cmt)
  
  out <- pkmod %>%
    param(
      ka = ka,
      CL = clin_CL, 
      V = clin_V, 
      PC = PC) %>%
    ev(regimen) %>%
    mrgsim(end = inter * ndoses, delta = 0.1)
  
  as_tibble(out)
}


plot_pk <- function(sim_df, MIC = NA, macIC90 = NA, casMBC90 = NA, ec50 = NA) {
  
  p <- ggplot(sim_df, aes(x = time, y = C_CP)) +
    geom_line() +
    labs(x = "Time (h)", y = "Concentration (mg/L)") +
    scale_y_log10() +
    theme_bw()
  
  if (!is.na(MIC)) {
    p <- p +
      geom_hline(yintercept = MIC, color = "red", linetype = "dashed") +
      geom_text(aes(x = Inf, y = MIC), 
                label = "MIC", 
                vjust = -1, 
                hjust = 1.1,
                color = "red")
  }
  
  if (length(macIC90) > 0 && !is.na(macIC90)) {
    p <- p +
      geom_hline(yintercept = macIC90, color = "blue", linetype = "dashed") +
      geom_text(aes(x = 1, y = macIC90), label = "macIC90", vjust = -1, color = "blue")
  }
  
  if (!is.na(casMBC90)) {
    p <- p +
      geom_hline(yintercept = casMBC90, color = "orange", linetype = "dashed") +
      geom_text(aes(x = 1, y = casMBC90), label = "casMBC90", vjust = -1, color = "orange")
  }
  
  if (!is.na(ec50)) {
    p <- p +
      geom_hline(yintercept = ec50, color = "purple", linetype = "dashed") +
      geom_text(aes(x = 1, y = ec50), label = "EC50", vjust = -1, color = "purple")
  }
  
  p
}


plot_lesion_pk <- function(sim_df, MIC = NA, macIC90 = NA, casMBC90 = NA, ec50 = NA) {
  
  p <- ggplot(sim_df, aes(x = time, y = C_LESION)) +
    geom_line() +
    labs(x = "Time (h)", y = "Concentration (mg/L)") +
    scale_y_log10() +
    theme_bw()
  
  if (!is.na(casMBC90)) {
    p <- p +
      geom_hline(yintercept = casMBC90, color = "orange", linetype = "dashed") +
      geom_text(aes(x = Inf, y = casMBC90), 
                label = "casMBC90", 
                vjust = 2, 
                hjust = 1.1,
                color = "orange")
  }
  
  if (!is.na(ec50)) {
    p <- p +
      geom_hline(yintercept = ec50, color = "purple", linetype = "dashed") +
      geom_text(aes(x = Inf, y = ec50), 
                label = "EC50", 
                vjust = 2, 
                hjust = 1.1, 
                color = "purple")+
      geom_hline(yintercept = ec50*9, color = "blue", linetype = "dashed") +
      geom_text(aes(x = Inf, y = ec50*9), 
                label = "EC90", 
                vjust = 2, 
                hjust = 1.1, 
                color = "blue")
    }
  
  p
}




# Run coverage simulation
run_coverage_sim <- function(clin_CL, clin_V, dose_vec, ndoses, inter,
                             MIC = NA, macIC90 = NA, casMBC90 = NA, ec50 = NA) {
  
  evs <- lapply(seq_along(dose_vec), function(i) {
    ev(
      ID   = i,
      amt  = dose_vec[i],
      ii   = inter,
      addl = ndoses - 1,
      cmt  = "CP"
    )
  })
  
  ev_all <- do.call(c, evs)
  
  out <- pkmod %>%
    param(CL = clin_CL, 
          V = clin_V) %>%
    ev(ev_all) %>%
    mrgsim(end = inter * ndoses, delta = 0.5)
  
  out_df <- as_tibble(out)
  
  out_df <- out_df %>%
    filter(time > 0) %>%
    mutate(dose = dose_vec[ID])   

  out_df %>%
    dplyr::filter(time %in% seq(24 * (ndoses - 1), 24 * ndoses, 0.1)) %>%
    dplyr::mutate(
      TAM = if (!is.null(MIC) && length(MIC) == 1 && !is.na(MIC)) {
        as.integer(C_CP > MIC)
      } else {
        0L
      },
      
      TAMAC = if (!is.null(macIC90) && length(macIC90) == 1 && !is.na(macIC90)) {
        as.integer(C_CP > macIC90)
      } else {
        0L
      },
      
      TACAS = if (!is.null(casMBC90) && length(casMBC90) == 1 && !is.na(casMBC90)) {
        as.integer(C_CP > casMBC90)
      } else {
        0L
      },
      
      TAE = if (!is.null(ec50) && length(ec50) == 1 && !is.na(ec50)) {
        as.integer(C_CP > ec50)
      } else {
        0L
      }
    ) %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(
      perc_TAM   = mean(TAM)   * 100,
      perc_TAMAC = mean(TAMAC) * 100,
      perc_TACAS = mean(TACAS) * 100,
      perc_TAE   = mean(TAE)   * 100
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(dose, .keep_all = TRUE)
}



run_coverage_sim_lesion <- function(clin_CL, clin_V, dose_vec, ndoses, inter, PC, KPL = 10,
                             MIC = NA, macIC90 = NA, casMBC90 = NA, ec50 = NA) {
  
  evs <- lapply(seq_along(dose_vec), function(i) {
    ev(
      ID   = i,
      amt  = dose_vec[i],
      ii   = inter,
      addl = ndoses - 1,
      cmt  = "CP"
    )
  })
  
  ev_all <- do.call(c, evs)
  
  out <- pkmod %>%
    param(CL = clin_CL, 
          V = clin_V,
          PC = PC) %>%
    ev(ev_all) %>%
    mrgsim(end = inter * ndoses, delta = 0.5)
  
  out_df <- as_tibble(out)
  
  out_df <- out_df %>%
    filter(time > 0) %>%
    mutate(dose = dose_vec[ID])   
  
  out_df %>%
    dplyr::filter(time %in% seq(24 * (ndoses - 1), 24 * ndoses, 0.1)) %>%
    dplyr::mutate(
      TAM = if (!is.null(MIC) && length(MIC) == 1 && !is.na(MIC)) {
        as.integer(C_LESION > MIC)
      } else {
        0L
      },
      
      TACAS = if (!is.null(casMBC90) && length(casMBC90) == 1 && !is.na(casMBC90)) {
        as.integer(C_LESION > casMBC90)
      } else {
        0L
      },
      
      TAE = if (!is.null(ec50) && length(ec50) == 1 && !is.na(ec50)) {
        as.integer(C_LESION > ec50)
      } else {
        0L
      },
      
      TAEC90 = if (!is.null(ec50) && length(ec50) == 1 && !is.na(ec50)) {
        as.integer(C_LESION > ec50*9)
      } else {
        0L
      }
      
      
    ) %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(
      perc_TAM    = mean(TAM)   * 100,
      perc_TACAS  = mean(TACAS) * 100,
      perc_TAE    = mean(TAE)   * 100,
      perc_TAEC90 = mean(TAEC90)   * 100
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(dose, .keep_all = TRUE)
  }


plot_coverage <- function(coverage_df, MIC = NA, macIC90 = NA, casMBC90 = NA, ec50 = NA) {
  p <- ggplot(coverage_df) +
    #scale_x_log10() +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100)) +
    labs(x = "Dose (mg)", y = "Coverage above target (%)") +
    theme_bw() +
    theme(legend.position = "none")
  
  if (!is.na(MIC)) {
    p <- p +
      geom_line(aes(x = dose, y = perc_TAM), color = "red", linewidth = 0.8, alpha = 0.8) +
      geom_text(
        data = coverage_df %>% filter(dose == max(dose)),
        aes(x = dose + 0.05 * (max(dose) - min(dose)), y = perc_TAM, label = "MIC"),
        vjust = -1, color = "red", show.legend = FALSE
      )
  }
  
  if (!is.na(casMBC90)) {
    p <- p +
      geom_line(aes(x = dose, y = perc_TACAS), color = "orange", linewidth = 0.8, alpha = 0.8) +
      geom_text(
        data = coverage_df %>% filter(dose == max(dose)),
        aes(x = dose + 0.05 * (max(dose) - min(dose)), y = perc_TACAS, label = "casMBC90"),
        vjust = -1, 
        hjust = 1,
        color = "orange", 
        show.legend = FALSE
      )
    }
  
  if (!is.na(ec50)) {
    p <- p +
      geom_line(aes(x = dose, y = perc_TAE), color = "purple", linewidth = 0.8, alpha = 0.8) +
      geom_text(
        data = coverage_df %>% filter(dose == max(dose)),
        aes(x = dose + 0.05 * (max(dose) - min(dose)), y = perc_TAE, label = "EC50"),
        vjust = -1, 
        hjust = 1,
        color = "purple", 
        show.legend = FALSE
      ) +
      geom_line(aes(x = dose, y = perc_TAEC90), color = "blue", linewidth = 0.8, alpha = 0.8) +
      geom_text(
        data = coverage_df %>% filter(dose == max(dose)),
        aes(x = dose + 0.05 * (max(dose) - min(dose)), y = perc_TAEC90, label = "EC90"),
        vjust = -1, 
        hjust = 1,
        color = "blue", 
        show.legend = FALSE
      )
  }
  
  p
}

calc_ec50_interp <- function(coverage_df, target = 50, coverage_col) {
  # Select doses and coverage percentages
  doses <- coverage_df$dose
  coverage <- coverage_df[[coverage_col]]
  
  if (all(is.na(coverage))) {
    return("no_input")
  }
  
  # Check if target is within range of coverage
  if (all(coverage < target)) {
    return(">10000")  # target never reached
  }
  
  if (all(coverage > target)) {
    return(min(doses))  # coverage always above target, return lowest dose
  }
  
  # Find indices where coverage crosses the target
  idx <- which(diff(coverage >= target) != 0)
  
  # If no crossing found (shouldn't happen due to above checks)
  if(length(idx) == 0) return("-")
  
  # We'll interpolate between doses[idx] and doses[idx+1]
  x1 <- doses[idx]
  x2 <- doses[idx + 1]
  y1 <- coverage[idx]
  y2 <- coverage[idx + 1]
  
  # Linear interpolation formula
  slope <- (y2 - y1) / (x2 - x1)
  dose_interp <- x1 + (target - y1) / slope
  
  return(as.integer(round(dose_interp)))
}


# Clinical EBA
code <- '
  $PARAM
    clin_CL   = 1,
    clin_Vss  = 10,
    EMAX      = 1,
    EC50      = 1,
    GAM       = 1,
    KIN       = 0,
    BACKBONE  = 0,
    base_log  = 6
  
  $CMT
    CENT      // drug amount
    BAC       // bacterial burden
    AUC       // cumulative AUC
    BACINT    // cumulative bacterial burden
  
  $OMEGA @labels nCL nVss nBAC
    0.09 0.09 
    0.09
  
  $MAIN
    double CL_i   = clin_CL * exp(nCL);
    double Vss_i  = clin_Vss * exp(nVss);
    
    if (NEWIND <= 1) {
      double log_start = base_log + nBAC; 
      BAC_0 = pow(10, log_start); 
    }
    
  $ODE
    // Plasma concentration
    double CP = CENT / Vss_i;
  
    // Emax model
    double EFF = BACKBONE + (EMAX * pow(CP, GAM)) / (pow(EC50, GAM) + pow(CP, GAM));
                 
    double log10CFU = log10(fmax(BAC, 1.0));
  
    // PK
    dxdt_CENT = -(CL_i / Vss_i) * CENT;
  
    // Bacterial dynamics
    dxdt_BAC = KIN * BAC - EFF * BAC;
  
    // Integrals
    dxdt_AUC    = CP;
    dxdt_BACINT = BAC;
  
  $CAPTURE
    CP       = CP;
    EFF      = EFF;
    log10CFU = log10CFU;
  '

eba_mod <- mcode("clinical_eba", code)


run_eba_sim <- function (clin_CL, clin_V, dose, ndoses, inter, EC50, Emax, GAM, backbone){
  
  regimen <- ev(
    amt = dose, 
    ii = inter, 
    addl = ndoses-1,
    cmt  = "CENT")
  
  pop_data <- data.frame(ID = 1:100)
  
  init_vals <- c(CENT = 0, BAC = 1e06, AUC = 0, BACINT = 0)
  
  out <- eba_mod %>%
    param(clin_CL = clin_CL, clin_Vss = clin_V, EMAX = Emax, EC50 = EC50, GAM = GAM, BACKBONE = backbone) %>%
    idata_set(pop_data) %>%
    ev(regimen) %>%
    init(BAC = 1e06) %>%   
    mrgsim(end = inter * ndoses, delta = 4) %>% # Smaller delta for smoother ribbons
    as.data.frame()
  
  #print(out)
  return(out)
  } 



plot_eba <- function(sim_df, drug_name, ndoses, inter){

  # Calculate Statistics
  summary_df <- sim_df %>%
    mutate(t_days = time / 24) %>%
    group_by(t_days) %>%
    summarise(
      median_BAC = median(BAC),
      low_BAC    = quantile(BAC, 0.05),
      high_BAC   = quantile(BAC, 0.95)
    )
  
  # 2. Reference Line Data
  x_end <- (inter * ndoses) / 24
  ref_label <- "RIF 600 mg"
  
  # Create a data frame for the reference line to allow mapping
  ref_df <- data.frame(
    x = 0, y = 10^6, 
    xend = (inter * ndoses) / 24, 
    yend = 10^(6 + (-0.176 * (inter * ndoses) / 24)),
    label = "RIF 600 mg"
  )
  
  sim_df <- sim_df %>%
    mutate(
      t_days = time / 24,
      BAC = pmax(BAC, 1e-8))  # avoid log10(0)
  
  ggplot(summary_df, aes(x = t_days)) +
    # Shaded Ribbon for IIV
    geom_ribbon(aes(ymin = low_BAC, ymax = high_BAC, fill = drug_name), alpha = 0.2) +
    # Median Line
    geom_line(aes(y = median_BAC, color = drug_name), linewidth = 1) +
    
    # Reference Line (mapped to the same color aesthetic to join the legend)
    geom_segment(data = ref_df, aes(x = x, y = y, xend = xend, yend = yend, linetype = label), 
                 color = "black", linewidth = 0.8) +
    
    scale_color_manual(values = setNames(c("red", "black"), c(drug_name, ref_label))) +
    #scale_fill_manual(values = setNames(c("red", "transparent"), c(drug_name, ref_label))) +
    
    labs(y = "Log10(CFU/mL)", 
         x = "Time (days)") +
    scale_x_continuous(breaks = seq(0, (inter*ndoses)/24, inter/24)) +
    scale_y_log10(breaks = 10^(0:8), labels = 0:8, limits = c(1, 1e8)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank()
          )
          
}


run_tcc_sim <- function(clin_CL, clin_Vss, dose, ndoses, inter, EC50, Emax, GAM, backbone, n_ind = 50, delta = 4) {
  
  regimen <- ev(
    amt  = dose,
    ii   = inter,
    addl = ndoses - 1,
    cmt  = "CENT"
  )
  
  init_vals <- c(CENT = 0, 
                 BAC = 1e6, 
                 AUC = 0, 
                 BACINT = 0) ### need to change initial value
  
  out <- eba_mod %>%
    param(
      clin_CL  = clin_CL,
      clin_Vss = clin_Vss,
      EMAX     = Emax,
      EC50     = EC50,
      GAM      = GAM,
      BACKBONE = backbone
    ) %>%
    ev(regimen) %>%
    init(init_vals) %>%
    mrgsim(nid = n_ind, 
           end = 56*24, 
           delta = delta)
  
  #print(names(as.data.frame(out)))
  #print(out)
  return(out)
}



plot_tcc <- function(sim_df, drug_name, max_days = 56) {
  
  print(sim_df)
  # --- 1. REFERENCE DATA PREP ---
  zenix_fit <- survfit(Surv(TIME, EVENT) ~ 1, data = bpal_zenix_tte)
  
  clinical_plot_data <- tidy(zenix_fit) %>%
    add_row(time = 0, estimate = 1, conf.low = 1, conf.high = 1, .before = 1) %>%
    mutate(time_end = lead(time)) %>%
    filter(!is.na(time_end))
  
  # --- 2. SIMULATION DATA PREP ---
  sim_km_df <- sim_df %>%
    as.data.frame() %>%
    mutate(t_days = .data$time / 24) %>%
    group_by(ID) %>%
    summarize(
      event_day = ifelse(any(BAC <= 1), min(t_days[BAC <= 1]), max_days),
      event = as.numeric(any(BAC <= 1)),
      .groups = "drop"
    ) 
  
  #print(head(sim_km_df))
  
  # Fit simulation KM curve
  sim_fit <- survfit(Surv(event_day, event) ~ 1, data = sim_km_df)
  sim_fit$call$data <- sim_df
  
  # 1. Transform the survfit object into a data frame
  # We use 'surv' (estimate), 'lower', and 'upper'
  raw_plot_data <- data.frame(
    time  = c(0, sim_fit$time),
    surv  = c(1, sim_fit$surv),
    lower = c(1, sim_fit$lower),
    upper = c(1, sim_fit$upper)
  )
  
  # 2. Create the "Step" data manually
  # This duplicates rows to create the 'corners' of the steps
  step_data <- raw_plot_data %>%
    arrange(time) %>%
    mutate(
      # Create the horizontal 'lag' values
      time_end = lead(time),
      surv_end = surv,
      lower_end = lower,
      upper_end = upper
    ) %>%
    # Filter out the last row which has no 'lead'
    filter(!is.na(time_end))
  
  ref_label <- "Zenix (BPaL)"
  
  ggplot() +
    # --- CLINICAL LAYER (Zenix) ---
    geom_rect(data = clinical_plot_data,
              aes(xmin = time, xmax = time_end, ymin = conf.low, ymax = conf.high, fill = ref_label),
              alpha = 0.1) +
    geom_step(data = clinical_plot_data,
              aes(x = time, y = estimate, color = ref_label),
              linewidth = 0.8) +
    
    # --- SIMULATION LAYER ---
    geom_rect(data = step_data, 
              aes(xmin = time, xmax = time_end, ymin = lower, ymax = upper, fill = drug_name),
              alpha = 0.2) +
    geom_step(data = step_data,
              aes(x = time, y = surv, color = drug_name),
              linewidth = 1) +
    
    # --- FORMATTING ---
    scale_color_manual(
      values = setNames(c("black", "#E41A1C"), c(ref_label, drug_name))
      ) +
    # 2. Fill scale for ribbons
    scale_fill_manual(
      values = setNames(c("grey20", "#E41A1C"), c(ref_label, drug_name))
    ) +
    scale_x_continuous(breaks = seq(0, 56, by = 7), 
                       limits = c(0, 56), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = "Time (Days)",
         y = "Proportion without Culture Conversion") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank())
}