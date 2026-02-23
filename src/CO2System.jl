module CO2System

"""
Julia adaptation of BFM CarbonateSystem routine based on MOCSY 2.0.
Module providing function CarbonateSystem to solve carbonate system variables using total alkalinity and
dissolved inorganic carbon concentrations
Uses SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use local Ks instead of its own

## INPUT:
  T    = in situ temperature                  [degrees C]
  S    = practical salinity                   [unitless]
  TA      = total alkalinity                     [ueq/kg]
  tc      = dissolved inorganic carbon           [umol/kg]
  pt      = total dissolved inorganic phosphorus [mmol/m3]
  sit     = total dissolved inorganic silicon    [mmol/m3]
  Bt      = total dissolved inorganic boron      computed
  St      = total dissolved inorganic sulfur     computed
  Ft      = total dissolved inorganic fluorine   computed
  K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  Patm    = atmospheric pressure [mbar=hPa]
  Rho     = in-situ densiaty     [kg/m3]
  p0   = hydrostatic pressure [dbar]
  If p0 is given, considers in situ T & total pressure (atm + hydrostatic) to compute fCO2 and pCO2
  otherwise, in situ T & only atm pressure (hydrostatic=0) for 'zero order' fCO2 and pCO2.

  ---------

## OUTPUT:
  ph   = pH on total scale
  pco2 = CO2 partial pressure (uatm)
  fco2 = CO2 fugacity (uatm)
  co2  = aqueous CO2 concentration in [umol/kg]
  hco3 = bicarbonate (HCO3-) concentration in [umol/kg]
  co3  = carbonate (CO3--) concentration in [umol/kg]
  OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  OmegaC = Omega for calcite, i.e., the   calcite saturation state

## ORIGINAL REFERENCES
  Orr and Epitaloni, 2015: "Improved routines to model the ocean
  carbonate system: mocsy 2.0." GMD 8.3 : 485-499.

"""

FT = Float64

mutable struct CarbonateSystemEquilibriumConstants
  K0 :: FT # solubility : [Co2]=k0Ac Pco2
  K1 :: FT # carbonate equilibrium I
  K2 :: FT # carbonate equilibrium II
  Kw :: FT # water dissociation
  Kb :: FT # constant for Boron equilibrium
  Ks :: FT # constant for bisulphate equilibrium
  Kf :: FT # constant for hidrogen fluoride equilibirum
  K1p :: FT
  K2p :: FT
  K3p :: FT # constants for phosphate equilibirum
  Ksi :: FT # constant for silicic acid equilibrium
  Kspc :: FT # constant for calcite equilibrium
  Kspa :: FT # constant for aragonite equilibrium
end
CarbonateSystemEquilibriumConstants() = CarbonateSystemEquilibriumConstants(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

const MEG :: FT = 1.e6
const PERMIL :: FT =1.e-3
const PERMEG :: FT = 1.e-6
const ZERO_KELVIN :: FT = -273.15
const MW_C :: FT = 12.011
const MW_Ca :: FT = 40.078
const Rgas :: FT = 83.14472
const Rgas_atm :: FT = 82.05736    # (cm3 * atm) / (mol * K)  CODATA (2006)
const vbarCO2 :: FT = 32.3        # partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
const CaRelCon :: FT = 0.02128     # Calcium ion relative concentration See Dickson (2007) and Munhoven (2013)
const p_atm0 :: FT = 1013.25 # Atmospheric pressure at sea level [mbar]
const bar2atm :: FT = 1000.0/p_atm0 # Conversion factor from bar to atm

export CarbonateSystem

"""
Function tosolve carbonate system variables using total alkalinity and
dissolved inorganic carbon concentrations

Arguments:
    S::FT: Salinity [psu]
    T::FT: Temperature [degrees C]
    Rho::FT: Density [kg/m^3]
    PO4::FT: Nitrate concentration [mmol/m3]
    Sil::FT: Silicate concentration [mmol/m3]
    DIC::FT: Dissolved inorganic carbon concentration [umol/kg]
    TA::FT: Total alkalinity [umol/kg]

Optional Arguments:
    patm::FT: Atmospheric pressure [mbar], default is p_atm0
    p0::FT: Pressure increment [mbar], default is 0.
    pH0::FT: Initial pH guess [pH], default is -1., so it will be calculated if not provided
    maxit::Int: Maximum number of iterations, default is 50

Returns:
    A tuple containing the calculated pH, CO2[umol/kg], HCO3[umol/kg], CO3[umol/kg], pCO2[uatm], OmegaC, OmegaA, fCO2[uatm], and an errorflag for the solution of the carbonate system equations (1: error, 0: no error).
"""
function CarbonateSystem(; S::FT,T::FT,Rho::FT,PO4::FT,Sil::FT,DIC::FT,TA::FT,patm::FT=p_atm0,p0::FT=0.0, pH0::FT=-1.0, maxit=50)

  CSeq = CarbonateSystemEquilibriumConstants()

  a0 =  (FT)[25.50, 15.82, 29.48, 25.60, 18.03,  9.78,
            48.76, 46.00, 14.51, 23.12, 26.57, 29.48]
  a1 =  (FT)[0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090,
            0.5304,  0.5304, 0.1211, 0.1758, 0.2020,  0.1622]
  a2 =  (FT)[0.0, 0.0, 2.608, -1.409, 0.316, -0.942,
            0.0, 0.0,-0.321, -2.647, -3.042, -2.608]
  b0 =  (FT)[3.08, -1.13, 2.84, 5.13, 4.53, 3.91,
            11.76, 11.76, 2.67, 5.15, 4.08, 2.84]
  b1 =  (FT)[0.0877, -0.1475, 0.0,  0.0794, 0.09,  0.054,
            0.3692,  0.3692, 0.0427, 0.09, 0.0714, 0.0  ]
  # b2 =  zeros(FT,12) unused

  errorflag = 0

  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # 1. PREPARE INPUT FIELDS
  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Convert input fields
  # from umol/kg to mol/kg
  ta = TA * PERMEG
  tc = DIC * PERMEG
  # from mmol/m^3 -> mol/kg
  pt  = PO4 / Rho * PERMIL
  sit = Sil / Rho * PERMIL

  # Absolute temperature (Kelvin) and related values
  tk     = T - ZERO_KELVIN  # T is in degC; tk is in degK
  temp2  = T^2
  tk100  = tk / 100.
  tk1002 = tk100^2
  invtk  = 1. / tk
  dlogtk = log(tk)

  # Salinity and simply related values
  s  = S
  s2 = S^2
  sqrts = sqrt(S)
  s15 = S^1.5

  # Hydrostatic Pressure [dbar]
  press = p0 * 0.1  # convert from dbar to bar
  pr2   = press^2 / Rgas
  pr    = press / Rgas

  # Atmospheric pressure
  pratm = patm * PERMIL  # convert from mbar to bar
  Ptot  = (pratm + press) * bar2atm   # total pressure (in atm) = atmospheric pressure + hydrostatic pressure

  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # 2. COMPUTE CONSTANTS
  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # chlorinity
  scl = s/1.80655

  # ionic strength
  is = 19.924*s/(1000.0-1.005*s)
  is2 = is^2
  sqrtis = sqrt(is)

  # Total concentrations for sulfate, fluoride, and boron
  # Sulfate: Morris & Riley (1966)
  St = 0.14 * scl/96.062
  # Fluoride: Riley (1965)
  Ft = 0.000067 * scl/18.9984
  # Boron: Lee et al (2010)
  Bt = 0.0002414 * scl/10.811

  # -----------------------------------------------------------------------
  # K0, solubility of co2 in the water (K Henry) from Weiss 1974
  # K0 = [CO2]/ fCO2 [mol/kg/atm]
  # -----------------------------------------------------------------------
  lnK = 93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
       s * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002)
  CSeq.K0 = exp( lnK )

  # -----------------------------------------------------------------------
  # Choice of Acidity constants
  # K1 = [H][HCO3]/[H2CO3]   ,   K2 = [H][CO3]/[HCO3]
  # -----------------------------------------------------------------------
  # Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
  CSeq.K1 = 10.0^(-1.0*(3633.86*invtk - 61.2172 + 9.6777*dlogtk -
       0.011555 * s + 0.0001152 * s2))
  CSeq.K2 = 10.0^(-1.0*(471.78*invtk + 25.9290 -
       3.16967*dlogtk - 0.01781 * s + 0.0001122 * s2))
  # Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
  # K1 = 10.0^(-1.0*( (6320.813*invtk + 19.568224*dlogtk -126.34048 +
  #      13.4038*sqrts + 0.03206*s - (5.242e-5)*s2) +
  #      (-530.659*sqrts - 5.8210*s)*invtk -2.0664*sqrts*dlogtk) )
  # K2 = 10.0^(-1.0*( (5143.692*invtk + 14.613358*dlogtk -90.18333 +
  #      21.3728*sqrts + 0.1218*s - (3.688e-4)*s2 ) +
  #      (-788.289*sqrts - 19.189*s)*invtk -3.374*sqrts*dlogtk) )
  # Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
  # K1 = 10.0^(-1.0*( (6320.813*invtk + 19.568224*dlogtk -126.34048 +
  #      13.409160*sqrts + 0.031646*s - (5.1895e-5)*s2 ) +
  #      (-531.3642*sqrts - 5.713*s)*invtk -2.0669166*sqrts*dlogtk) )
  # K2 = 10.0^(-1.0*( (5143.692*invtk + 14.613358*dlogtk -90.18333 +
  #      21.225890*sqrts + 0.12450870*s - (3.7243e-4)*s2 ) +
  #      (-779.3444*sqrts - 19.91739*s)*invtk -3.3534679*sqrts*dlogtk) )

  #-----------------------------------------------------------------------
  # Kb = [H][BO2]/[HBO2]
  # Millero p.669 (1995) using data from Dickson (1990)    (total scale)
  #-----------------------------------------------------------------------
  lnK = (-8966.90 - 2890.53*sqrts - 77.942*s +
       1.728*s15 - 0.0996*s2)*invtk +
       (148.0248 + 137.1942*sqrts + 1.62142*s) +
       (-24.4344 - 25.085*sqrts - 0.2474*s) *
       dlogtk + 0.053105*sqrts*tk
  CSeq.Kb = exp(lnK)

  # -----------------------------------------------------------------------
  # K1p = [H][H2PO4]/[H3PO4]
  # DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  # Millero (1995), p.670, eq. 65                        (seawater scale)
  # -----------------------------------------------------------------------
  lnK = -4576.752*invtk + 115.540 - 18.453 * dlogtk +
       (-106.736*invtk + 0.69171) * sqrts +
       (-0.65643*invtk - 0.01844) * s
  CSeq.K1p = exp(lnK)

  # -----------------------------------------------------------------------
  # K2p = [H][HPO4]/[H2PO4]
  # DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  # Millero (1995), p.670, eq. 66                        (seawater scale)
  # -----------------------------------------------------------------------
  lnK = -8814.715*invtk + 172.1033 - 27.927 * dlogtk +
       (-160.340*invtk + 1.3566) * sqrts +
       (0.37335*invtk - 0.05778) * s
  CSeq.K2p = exp(lnK)

  #------------------------------------------------------------------------
  # K3p = [H][PO4]/[HPO4]
  # DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  # Millero (1995), p.670, eq. 67                        (seawater scale)
  # -----------------------------------------------------------------------
  lnK = -3070.75*invtk - 18.126 +
       (17.27039*invtk + 2.81197) *
       sqrts + (-44.99486*invtk - 0.09984) * s
  CSeq.K3p = exp(lnK)

  #------------------------------------------------------------------------
  # ksi = [H][SiO(OH)3]/[Si(OH)4]
  # Millero (1995), p.671, eq. 72                        (seawater scale)
  # -----------------------------------------------------------------------
  lnK = -8904.2*invtk + 117.400 - 19.334 * dlogtk +
       (-458.79*invtk + 3.5913) * sqrtis +
       (188.74*invtk - 1.5998) * is +
       (-12.1652*invtk + 0.07871) * is2 +
       log(1.0-0.001005*s)
  CSeq.Ksi = exp(lnK)

  # -----------------------------------------------------------------------
  # Kw = [H][OH]
  # Millero (1995) p.670, eq. 63 from composite data     (seawater scale)
  # -----------------------------------------------------------------------
  lnK = 148.9802 -13847.26*invtk - 23.6521 * dlogtk +
       (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
       sqrts - 0.01615 * s
  CSeq.Kw = exp(lnK)

  #------------------------------------------------------------------------
  # ks = [H][SO4]/[HSO4]
  # Dickson (1990, J. chem. Thermodynamics 22, 113)          (free scale)
  #------------------------------------------------------------------------
  lnK = -4276.1*invtk + 141.328 - 23.093*dlogtk +
       (-13856.0*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
       (35474.0*invtk - 771.54 + 114.723*dlogtk) * is -
       2698.0*invtk*is^1.5 + 1776.0*invtk*is2 +
       log(1. - 0.001005*s)
  Ks_0p = exp(lnK)

  #------------------------------------------------------------------------
  # kf = [H][F]/[HF]
  # Perez & Fraga (1987) recom. by Dickson et al., (2007)   (total scale)
  #------------------------------------------------------------------------
  lnK = 874.0*invtk - 9.68 + 0.111*sqrts
  Kf_0p = exp(lnK)

  #------------------------------------------------------------------------
  # Kspc = [Ca2+] [CO32-] - apparent solubility product of Calcite
  # Mucci (1983)  [mol/kg-soln]
  #------------------------------------------------------------------------
  CSeq.Kspc = 10.0^(1.0*( -171.9065 - 0.077993 * tk + 2839.319 * invtk +
          71.595 * log10(tk) + sqrts * (-0.77712 +
          0.0028426 * tk + 178.34 * invtk) -
          0.07711 * s + 0.0041249 *s15 ) )
  #------------------------------------------------------------------------
  # Kspa = [Ca2+] [CO32-] - apparent solubility product of Aragonite
  # Mucci (1983)  [mol/kg-soln]
  #------------------------------------------------------------------------
  CSeq.Kspa = 10.0^(1.0*(  -171.945 - 0.077993 * tk + 2903.293 * invtk +
          71.595 * log10(tk) + sqrts * (-0.068393 +
          0.0017276 * tk + 88.135 * invtk) -
        0.10018 * s + 0.0059415 * s15 ) )

  # Pressure effect on K0 based on Weiss (1974, equation 5)
  CSeq.K0 = CSeq.K0 * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk) )

  # Pressure effect on all other K's (based on Millero, (1995)
  # Index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
  #           K1p(9), K2p(10), K3p(11), Ksi(12)
  if press > 0.
    deltav  =  -a0 .+ T*a1 .+ 1.e-3*temp2*a2
    deltak  = 1.0e-3*(-b0 .+ T*b1)
    lnkpok0 = -invtk*pr*deltav .+ 0.5*invtk*pr2*deltak
  else
    lnkpok0 = zeros(FT,12)
  end

  # Pressure correction on Ks (Free scale)
  CSeq.Ks = Ks_0p*exp(lnkpok0[5])
  # Conversion factor total -> free scale
  total2free     = 1.0/(1.0 + St/CSeq.Ks)   # Kfree = Ktotal*total2free
  # Conversion factor total -> free scale at pressure zero
  total2free_0p  = 1.0/(1.0 + St/Ks_0p)   # Kfree = Ktotal*total2free

  # Pressure correction on Kf
  # Kf must be on FREE scale before correction
  Kf_0p = Kf_0p * total2free_0p   # Convert from Total to Free scale (pressure 0)
  CSeq.Kf    = Kf_0p * exp(lnkpok0[6]) # Pressure correction (on Free scale)
  CSeq.Kf    = CSeq.Kf / total2free         # Convert back from Free to Total scale

  # Convert between seawater and total hydrogen (pH) scales
  free2SWS  = 1. + St/CSeq.Ks + Ft/(CSeq.Kf*total2free)  # using Kf on free scale
  total2SWS = total2free * free2SWS                  # KSWS = Ktotal*total2SWS
  SWS2total = 1. / total2SWS
  # Conversion at pressure zero
  free2SWS_0p  = 1. + St/Ks_0p + Ft/(Kf_0p)  # using Kf on free scale
  total2SWS_0p = total2free_0p * free2SWS_0p         # KSWS = Ktotal*total2SWS

  # Convert from Total to Seawater scale before pressure correction
  # Must change to SEAWATER scale: K1, K2, Kb
  CSeq.K1 = CSeq.K1 * total2SWS_0p
  CSeq.K2 = CSeq.K2 * total2SWS_0p
  CSeq.Kb = CSeq.Kb * total2SWS_0p
  # Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw
  # Other contants (keep on another scale):
  #    - K0         (independent of pH scale, already pressure corrected)
  #    - Ks         (already on Free scale;   already pressure corrected)
  #    - Kf         (already on Total scale;  already pressure corrected)
  #    - Kspc, Kspa (independent of pH scale; pressure-corrected below)

  # Perform actual pressure correction (on seawater scale)
  CSeq.K1   = CSeq.K1   * exp(lnkpok0[1])
  CSeq.K2   = CSeq.K2   * exp(lnkpok0[2])
  CSeq.Kb   = CSeq.Kb   * exp(lnkpok0[3])
  CSeq.Kw   = CSeq.Kw   * exp(lnkpok0[4])
  CSeq.Kspc = CSeq.Kspc * exp(lnkpok0[7])
  CSeq.Kspa = CSeq.Kspa * exp(lnkpok0[8])
  CSeq.K1p  = CSeq.K1p  * exp(lnkpok0[9])
  CSeq.K2p  = CSeq.K2p  * exp(lnkpok0[10])
  CSeq.K3p  = CSeq.K3p  * exp(lnkpok0[11])
  CSeq.Ksi  = CSeq.Ksi  * exp(lnkpok0[12])

  # Convert back to original total scale:
  CSeq.K1  = CSeq.K1  * SWS2total
  CSeq.K2  = CSeq.K2  * SWS2total
  CSeq.K1p = CSeq.K1p * SWS2total
  CSeq.K2p = CSeq.K2p * SWS2total
  CSeq.K3p = CSeq.K3p * SWS2total
  CSeq.Kb  = CSeq.Kb  * SWS2total
  CSeq.Ksi = CSeq.Ksi * SWS2total
  CSeq.Kw  = CSeq.Kw  * SWS2total

  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # 3. COMPUTE CARBONATE SYSTEM
  #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # H+ concentration (mol/kg) at previous step
  if pH0 <= 0.0
     Hi = Hini_for_at(CSeq,ta,tc,Bt)
  else
     Hi = 10.0^(-pH0)
  end

  # Solve for H+ using above result as the initial H+ value (mol/kg)
  H = solve_at_general(CSeq,ta, tc, Bt, pt, sit, St, Ft, Hi, maxit)

  errorflag = H < 0.0 ? 1 : errorflag

  # Calculate pH  from H+ concentration (mol/kg)
  pH = -1. * log10( H )

  # Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0 + CSeq.Ks*(1.0 + St/CSeq.Ks)/H)
  HF   = 1.0 / (1.0 + CSeq.Kf/H)
  HSI  = 1.0 / (1.0 + H/CSeq.Ksi)
  HPO4 = (CSeq.K1p*CSeq.K2p*(H + 2.0*CSeq.K3p) - H^3)    /
         (H^3 + CSeq.K1p*H^2 + CSeq.K1p*CSeq.K2p*H + CSeq.K1p*CSeq.K2p*CSeq.K3p)
  ab = Bt/(1.0 + H/CSeq.Kb)
  aw = CSeq.Kw/H - H/(1.0 + St/CSeq.Ks)
  ac = ta + HSO4 - sit*HSI - ab - aw + Ft*HF - pt*HPO4

  # Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0 * tc - ac) / (2.0 + CSeq.K1 / H)
  cb = CSeq.K1 * cu / H
  cc = CSeq.K2 * cb / H

  # Determine Omega Calcite and Aragonite (see Munhoven (2013, GMD))
  Ca = (CaRelCon / MW_Ca) * s/1.80655
  OmegaA = (Ca*cc) / CSeq.Kspa
  OmegaC = (Ca*cc) / CSeq.Kspc

  # Determine CO2 fugacity [uatm]
  fCO2 = cu * 1.e6/CSeq.K0

  # Determine CO2 partial pressure from CO2 fugacity [uatm]
  # compute fugacity coefficient terms : B, Del, xc2
  B = -1636.75 + 12.0408*tk - 0.0327957*(tk*tk) + 0.0000316528*(tk*tk*tk)
  Del = 57.7 - 0.118*tk
  xCO2approx = fCO2 * 1.e-6

  if press > 0
     xCO2approx *= exp( ((1.0-Ptot)*32.3)/(82.05736*tk) )   # of K0 press. correction, see Weiss (1974, equation 5)
  end

  xc2 = (1. - xCO2approx)^2
  fugcoeff = exp( Ptot*(B + 2.0*xc2*Del)/(Rgas_atm*tk) )
  pco2 = fCO2 / fugcoeff

  # scale from mol/kg -----> umol/kg
  co2  = cu * MEG
  hco3 = cb * MEG
  co3  = cc * MEG

  return (CO2=co2, HCO3=hco3, CO3=co3, pCO2=pco2, OmegaC=OmegaC, OmegaA=OmegaA, fCO2=fCO2, pH=pH, Hplus=H, errorflag=errorflag)

  end

# ----------------

  function Hini_for_at(cq,p_alkcb, p_dictot, p_bortot)
  # Function returns the root for the 2nd order approximation of the
  # DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
  # around the local minimum, if it exists.

  # Returns * 1E-03 if p_alkcb <= 0
  #         * 1E-10 if p_alkcb >= 2*p_dictot + p_bortot
  #         * 1E-07 if 0 < p_alkcb < 2*p_dictot + p_bortot
  #                    and the 2nd order approximation does not have a solution

    if p_alkcb <= 0.
      p_hini = 1.e-3
    elseif (p_alkcb >= (2.0*p_dictot + p_bortot))
      p_hini = 1.e-10
    else
      zca = p_dictot/p_alkcb
      zba = p_bortot/p_alkcb

      # Coefficients of the cubic polynomial
      za2 = cq.Kb*(1.0 - zba) + cq.K1*(1.0-zca)
      za1 = cq.K1*cq.Kb*(1.0 - zba - zca) + cq.K1*cq.K2*(1.0 - (zca+zca))
      za0 = cq.K1*cq.K2*cq.Kb*(1. - zba - (zca+zca))
                                            # Taylor expansion around the minimum
      zd = za2*za2 - 3.0*za1              # Discriminant of the quadratic equation
                                            # for the minimum close to the root

      if zd > 0.                    # If the discriminant is positive
        zsqrtd = sqrt(zd)
        if za2 < 0.
          zhmin = (-za2 + zsqrtd)/3.
        else
          zhmin = -za1/(za2 + zsqrtd)
        end
        p_hini = zhmin + sqrt(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
      else
        p_hini = 1.e-7
      end
    end

    Hini_for_at = p_hini

  end

function solve_at_general(cseq, p_alktot, p_dictot, p_bortot,
                              p_po4tot, p_siltot,
                              p_so4tot, p_flutot,
                              p_h, MaxIterPHsolver)

    # Purpose: Compute [H+] ion concentration from sea-water ion concentrations,
    #          alkalinity, DIC, and equilibrium constants
    # Universal pH solver that converges from any given initial value,
    # determines upper an lower bounds for the solution if required

    # local constants
    #-----------------
    pz_exp_threshold :: FT = 1.0
    pp_rdel_ah_target :: FT = 1.E-8


    # TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
    aphscale = 1. + p_so4tot/cseq.Ks

    # lower and upper bounds of "non-water-selfionization"
    # contributions to total alkalinity (the infimum and the supremum), i.e
    # inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])
    zalknw_inf = -p_po4tot - p_so4tot - p_flutot
    zalknw_sup = p_dictot + p_dictot + p_bortot + p_po4tot + p_po4tot + p_siltot

    zdelta = (p_alktot-zalknw_inf)^2 + 4.0*cseq.Kw/aphscale

    if p_alktot >= zalknw_inf
       zh_min = 2.0*cseq.Kw /( p_alktot-zalknw_inf + sqrt(zdelta) )
    else
       zh_min = aphscale*(-(p_alktot-zalknw_inf) + sqrt(zdelta) ) / 2.
    end

    zdelta = (p_alktot-zalknw_sup)^2 + 4.0*cseq.Kw/aphscale

    if p_alktot <= zalknw_sup
       zh_max = aphscale*(-(p_alktot-zalknw_sup) + sqrt(zdelta) ) / 2.
    else
       zh_max = 2.0*cseq.Kw /( p_alktot-zalknw_sup + sqrt(zdelta) )
    end

    zh = max(min(zh_max, p_h), zh_min)
    niter_atgen        = 0                 # Reset counters of iterations
    zeqn_absmin        = floatmax(FT)

    while true
       if niter_atgen >= MaxIterPHsolver
          zh = -1.0
          break
       end

       zh_prev = zh

       # Compute total alkalinity from ion concentrations and equilibrium constants
       # H2CO3 - HCO3 - CO3 : n=2, m=0
       znumer     = 2.0*cseq.K1*cseq.K2 + zh*       cseq.K1
       zdenom     =         cseq.K1*cseq.K2 + zh*(      cseq.K1 + zh)
       zalk_dic   = p_dictot * (znumer/zdenom)
       zdnumer    = cseq.K1*cseq.K1*cseq.K2 + zh*(4.0*cseq.K1*cseq.K2 + zh* cseq.K1 )
       zdalk_dic  = -p_dictot * (zdnumer/zdenom^2)
       # B(OH)3 - B(OH)4 : n=1, m=0
       znumer     = cseq.Kb
       zdenom     = cseq.Kb + zh
       zalk_bor   = p_bortot * (znumer/zdenom)
       zdnumer    = cseq.Kb
       zdalk_bor  = -p_bortot * (zdnumer/zdenom^2)
       # H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
       znumer     = 3.0*cseq.K1p*cseq.K2p*cseq.K3p + zh*(2.0*cseq.K1p*cseq.K2p + zh* cseq.K1p)
       zdenom     =         cseq.K1p*cseq.K2p*cseq.K3p + zh*(        cseq.K1p*cseq.K2p + zh*(cseq.K1p + zh))
       zalk_po4   = p_po4tot * (znumer/zdenom - 1.) # Zero level of H3PO4 = 1
       zdnumer    = cseq.K1p*cseq.K2p*cseq.K1p*cseq.K2p*cseq.K3p +
           zh*(4.0*cseq.K1p*cseq.K1p*cseq.K2p*cseq.K3p +
           zh*(9.0*cseq.K1p*cseq.K2p*cseq.K3p + cseq.K1p*cseq.K1p*cseq.K2p +
           zh*(4.0*cseq.K1p*cseq.K2p +
           zh*cseq.K1p)))
       zdalk_po4  = -p_po4tot * (zdnumer/zdenom^2)
       # H4SiO4 - H3SiO4 : n=1, m=0
       znumer     = cseq.Ksi
       zdenom     = cseq.Ksi + zh
       zalk_sil   = p_siltot * (znumer/zdenom)
       zdnumer    = cseq.Ksi
       zdalk_sil  = -p_siltot * (zdnumer/zdenom^2)
       # HSO4 - SO4 : n=1, m=1
       znumer     = cseq.Ks
       zdenom     = cseq.Ks + zh
       zalk_so4   = p_so4tot * (znumer/zdenom - 1.)
       zdnumer    = cseq.Ks
       zdalk_so4  = -p_so4tot * (zdnumer/zdenom^2)
       # HF - F : n=1, m=1
       znumer     = cseq.Kf
       zdenom     = cseq.Kf + zh
       zalk_flu   = p_flutot * (znumer/zdenom - 1.)
       zdnumer    = cseq.Kf
       zdalk_flu  = -p_flutot * (zdnumer/zdenom^2)
       # H2O - OH
       zalk_wat   = cseq.Kw/zh - zh/aphscale
       zdalk_wat  = -cseq.Kw/zh^2 - 1.0/aphscale

       zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil +
            zalk_so4 + zalk_flu + zalk_wat - p_alktot

       zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil +
            zdalk_so4 + zdalk_flu + zdalk_wat


       # Adapt bracketing interval
       if zeqn > 0.0
          zh_min = zh_prev
       elseif zeqn < 0.0
          zh_max = zh_prev
       else
          # zh is the root; unlikely but, one never knows
          break
       end
       # Now determine the next iterate zh
       niter_atgen += 1

       if abs(zeqn) >= 0.5*zeqn_absmin
          zh = sqrt(zh_max * zh_min)
          zh_lnfactor = (zh - zh_prev)/zh_prev # Required to test convergence below
       else
          zh_lnfactor = -zeqn/(zdeqndh*zh_prev)
          if abs(zh_lnfactor) > pz_exp_threshold
             zh          = zh_prev*exp(zh_lnfactor)
          else
             zh_delta    = zh_lnfactor*zh_prev
             zh          = zh_prev + zh_delta
          end
          if zh < zh_min # if [H]_new < [H]_min
             zh                = sqrt(zh_prev * zh_min)
             zh_lnfactor       = (zh - zh_prev)/zh_prev # Required to test convergence below
          end
          if zh > zh_max  # if [H]_new > [H]_max
             zh                = sqrt(zh_prev * zh_max)
             zh_lnfactor       = (zh - zh_prev)/zh_prev # Required to test convergence below
          end
       end
       zeqn_absmin = min( abs(zeqn), zeqn_absmin)

       (abs(zh_lnfactor) < pp_rdel_ah_target) && break

    end

    solve_at_general = zh

end

end
