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
model_object <- rxode2({
  # pre-calculation
  V_central = V1 * exp(eta_V) * BW 
  k_elimination = (CL * exp(eta_CL) * CL_creatinine) / V_central
  
  # ODEs
  d/dt(m_central) = -(k_elimination + k12) * m_central + k21 * m_peripheral
  d/dt(m_peripheral) = k12 * m_central - k21 * m_peripheral
  
  # post-calculation
  C_central = (m_central / V_central) + add.err # is this the correct syntax? nlmixr2 syntax does not seem to work
})