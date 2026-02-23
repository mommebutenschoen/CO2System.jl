# CO2System

Julia adaptation of the CarbonateSystem routine from the Biogeochemical Flux Model (Fortran) based on MOCSY 2.0.
This Module provides a function CarbonateSystem to solve carbonate system variables using total alkalinity and
dissolved inorganic carbon concentrations
Uses SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use local Ks instead of its own

### Input:

   S:practical salinity [unitless]  
   T: in situ temperature [degrees C]  
   Rho: in-situ densiaty [kg/m3]
   PO4: dissolved inorganic phosphorus [mmol/m3]  
   Sil: total dissolved inorganic silicon [mmol/m3]  
   DIC: dissolved inorganic carbon [umol/kg]  
   TA: total alkalinity [ueq/kg]  
   patm: atmospheric pressure [mbar=hPa]  
   p0   = hydrostatic pressure [dbar]

### Parameters computed internally
   Bt: total dissolved inorganic boron    
   St      = total dissolved inorganic sulfur  
   Ft      = total dissolved inorganic fluorine  
   K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  
  If p0 is given, considers in situ T & total pressure (atm + hydrostatic) to compute fCO2 and pCO2
  otherwise, in situ T & only atm pressure (hydrostatic=0) for 'zero order' fCO2 and pCO2.

  ---------

### Output:
  ph   = pH on total scale
  pco2 = CO2 partial pressure (uatm)
  fco2 = CO2 fugacity (uatm)
  co2  = aqueous CO2 concentration in [umol/kg]
  hco3 = bicarbonate (HCO3-) concentration in [umol/kg]
  co3  = carbonate (CO3--) concentration in [umol/kg]
  OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  OmegaC = Omega for calcite, i.e., the   calcite saturation state

### ORIGINAL REFERENCES
  Orr and Epitaloni, 2015: "Improved routines to model the ocean
  carbonate system: mocsy 2.0." GMD 8.3 : 485-499.

[![Build Status](https://github.com/mommebutenschoen/CO2System.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mommebutenschoen/CO2System.jl/actions/workflows/CI.yml?query=branch%3Amain)
