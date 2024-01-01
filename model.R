# LIBRARIES----------
library(rxode2)

# SUPPORTING FUNCTIONS---------
# Cited creatinine clearance translation from steady-state plasma creatinine concentration
creatinineCLModel <- function(age, bw, C_creatinine){
  # creatinine clearance estimation by Cockroft and Gault (1976; cited by Aarons in the core source publication)
  (150 - age) * bw / C_creatinine
}

# INITIAL STATE -----------
initialState <- c(m_central = 0, m_peripheral = 0)

# MODEL OBJECT # --------------
ab_amn_1989_tob_aarons <- rxode2({
  # pre-calculation
  V_central = V_ratio * exp(eta_V_ratio) * BW 
  k_elimination = (CL_sigma * exp(eta_CL_sigma) * CL_creatinine) / V_central
  
  # ODEs
  d/dt(m_central) = -(k_elimination + k12) * m_central + k21 * m_peripheral
  d/dt(m_peripheral) = k12 * m_central - k21 * m_peripheral
  
  # post-calculation
  C_central = m_central / V_central
})