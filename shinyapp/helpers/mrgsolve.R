library(ggplot2)
library(mrgsolve)

# Lesion PK model (Oral Dosing)
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
