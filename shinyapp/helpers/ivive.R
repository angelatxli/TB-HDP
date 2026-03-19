
species_params <- data.frame(
  Species = c("Human", "Mouse", "Rat", "Dog"),
  BW = c(60, 0.02, 0.25, 10),
  LiverWeight = c(25.7, 80, 40, 25),  # liver weight (g/kg)
  MPPGL = c(40, 50, 45, 30), # microsomal protein per gram liver (mg/g liver)
  HPGL = c(120, 135, 120, 100), # million cells per gram liver
  Qh = c(20.7, 126, 77.4, 56.1) # hepatic blood flow (mL/min/kg)
)



well_stirred_model <- function(CL_LM_invivo, fu, Qh) {
  CLh <- (Qh * fu * CL_LM_invivo) / (Qh + fu * CL_LM_invivo)
  return(CLh)
}


lm_ivive <- function(CL_LM, fu) {
  
  # Step 1: Scale in vitro CLint to in vivo (mL/min/kg)
  CL_LM_invivo = CL_LM*40*25.7*10^-3
  
  # Step 2: Apply well-stirred model
  CLh_LM_result = well_stirred_model(CL_LM_invivo, fu, 20.7)
  
  # Step 3: Convert to L/h (multiply by BW in kg and convert to L)
  CLh_LM_result_L_h = CLh_LM_result * 55 * 60/1000
  
  return(CLh_LM_result_L_h)
}



hep_ivive <- function(CL_LM, fu) {
  
  # Step 1: Scale in vitro CLint to in vivo (mL/min/kg)
  CL_LM_invivo = CL_LM*120*25.7*10^-3
  
  # Step 2: Apply well-stirred model
  CLh_LM_result = well_stirred_model(CL_LM_invivo, fu, 20.7)
  
  # Step 3: Convert to L/h (multiply by BW in kg and convert to L)
  CLh_LM_result_L_h = CLh_LM_result * 55 * 60/1000
  
  return(CLh_LM_result_L_h)
}

