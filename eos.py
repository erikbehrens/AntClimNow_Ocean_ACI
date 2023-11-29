#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 07:44:24 2023

@author: behrense
"""

from numpy import zeros, ones, sqrt, abs, floor, int32, reshape, nan, asarray
from numpy import mean, nanmax, sum, isfinite, isnan, sin, pi
from scipy.io import loadmat
from copy import deepcopy
from scipy.interpolate import interp1d

def sigmai_dep(CT, SA, p):
  """
  Compute the in-situ (or potential) density from CONSERVATIVE Temperature, 
  ABSOLUTE Salinity and pressure fields using the TEOS-10 EOS
  
  p can be a field (computed with subroutine p_from_z say) to get in-situ 
  density or can be a number p = p_ref, then it is a potential density 
  referenced to p_ref
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    CT             = Conservative Temperature      t        deg celsius
    SA             = Absolute Salinity             s        g / kg
    p              = (reference) pressure          p        dbar
    
  Returns:
    sigmai_dep_out = (potential density            rho      kg / m^3
    
  """
  # ensures that SA is non-negative.
  SA = abs(SA)

  # deltaS = 24
  sfac = 0.0248826675584615                 # sfac   = 1/(40*(35.16504/35)).
  offset = 5.971840214030754e-1             # offset = deltaS*sfac.

  x2 = sfac * SA
  xs = sqrt(x2 + offset)
  ys = CT * 0.025
  z  = p * 1e-4

  v000 =  1.0769995862e-3
  v001 = -6.0799143809e-5
  v002 =  9.9856169219e-6
  v003 = -1.1309361437e-6
  v004 =  1.0531153080e-7
  v005 = -1.2647261286e-8
  v006 =  1.9613503930e-9
  v010 = -1.5649734675e-5
  v011 =  1.8505765429e-5
  v012 = -1.1736386731e-6
  v013 = -3.6527006553e-7
  v014 =  3.1454099902e-7
  v020 =  2.7762106484e-5
  v021 = -1.1716606853e-5
  v022 =  2.1305028740e-6
  v023 =  2.8695905159e-7
  v030 = -1.6521159259e-5
  v031 =  7.9279656173e-6
  v032 = -4.6132540037e-7
  v040 =  6.9111322702e-6
  v041 = -3.4102187482e-6
  v042 = -6.3352916514e-8
  v050 = -8.0539615540e-7
  v051 =  5.0736766814e-7
  v060 =  2.0543094268e-7
  v100 = -3.1038981976e-4
  v101 =  2.4262468747e-5
  v102 = -5.8484432984e-7
  v103 =  3.6310188515e-7
  v104 = -1.1147125423e-7
  v110 =  3.5009599764e-5
  v111 = -9.5677088156e-6
  v112 = -5.5699154557e-6
  v113 = -2.7295696237e-7
  v120 = -3.7435842344e-5
  v121 = -2.3678308361e-7
  v122 =  3.9137387080e-7
  v130 =  2.4141479483e-5
  v131 = -3.4558773655e-6
  v132 =  7.7618888092e-9
  v140 = -8.7595873154e-6
  v141 =  1.2956717783e-6
  v150 = -3.3052758900e-7
  v200 =  6.6928067038e-4
  v201 = -3.4792460974e-5
  v202 = -4.8122251597e-6
  v203 =  1.6746303780e-8
  v210 = -4.3592678561e-5
  v211 =  1.1100834765e-5
  v212 =  5.4620748834e-6
  v220 =  3.5907822760e-5
  v221 =  2.9283346295e-6
  v222 = -6.5731104067e-7
  v230 = -1.4353633048e-5
  v231 =  3.1655306078e-7
  v240 =  4.3703680598e-6
  v300 = -8.5047933937e-4
  v301 =  3.7470777305e-5
  v302 =  4.9263106998e-6
  v310 =  3.4532461828e-5
  v311 = -9.8447117844e-6
  v312 = -1.3544185627e-6
  v320 = -1.8698584187e-5
  v321 = -4.8826139200e-7
  v330 =  2.2863324556e-6
  v400 =  5.8086069943e-4
  v401 = -1.7322218612e-5
  v402 = -1.7811974727e-6
  v410 = -1.1959409788e-5
  v411 =  2.5909225260e-6
  v420 =  3.8595339244e-6
  v500 = -2.1092370507e-4
  v501 =  3.0927427253e-6
  v510 =  1.3864594581e-6
  v600 =  3.1932457305e-5

  v = v000 + ( 
        xs * (v100 + xs * (v200 + xs * (v300 + xs * (v400 + xs * (v500 
      + v600 * xs))))) + ys * (v010 + xs * (v110 + xs * (v210 + xs * (v310 
      + xs * (v410 + v510 * xs)))) + ys * (v020 + xs * (v120 + xs * (v220 
      + xs * (v320 + v420 * xs))) + ys * (v030 + xs * (v130 + xs * (v230 
      + v330 * xs)) + ys * (v040 + xs * (v140 + v240*xs) + ys * (v050 
      + v150 * xs + v060 * ys))))) + z * (v001 + xs * (v101 + xs * (v201 
      + xs * (v301 + xs * (v401 + v501 * xs)))) + ys * (v011 + xs * (v111
      + xs * (v211 + xs * (v311 + v411 * xs))) + ys * (v021 + xs * (v121 
      + xs * (v221 + v321 * xs)) + ys * (v031 + xs * (v131 + v231 * xs) 
      + ys * (v041 + v141 * xs + v051 * ys)))) + z * (v002 + xs * (v102 
      + xs * (v202 + xs * (v302 + v402 * xs))) + ys * (v012 + xs * (v112 
      + xs * (v212 + v312 * xs)) + ys * (v022 + xs * (v122 + v222 * xs) 
      + ys * (v032 + v132 * xs + v042 * ys))) + z * (v003 + xs * (v103 
      + v203 * xs) + ys * (v013 + v113 * xs + v023 * ys) + z * (v004 
      + v104 * xs + v014 * ys + z * (v005 + v006 * z)))))
              )

  sigmai_dep_out = (1 / v) - 1000
  
  return sigmai_dep_out



def sigmantr(CT, SA, z, lon, lat):
  """
  Compute the neutral volumic mass (kg/m3) from known CONSERVATIVE Temperature 
  and PRACTICAL salinity fields from the Jackett and McDougall (2005) EOS
  
  Done via a conversion of CT to pt and SA to ps routines. The CT and SA fields
  should be massaged into a 3d field in (depth, lat, lon) arrangement even if
  one of the dimensions is only of length 1 so that the cycling through the
  depth index is ok
  
  Inputs:
    CT           = Conservative Temperature      t       deg celsius
    SA           = Absolute Salinity             s         g / kg
    z            = z-location of data                      m
                   (negative z is ocean depth)
    lon          = longitude of data                     deg
    lat          = latitude  of data                     deg
    
  Returns:
    sigmantr_out = neutral density               rho      kg / m^3
    
  """
  zrau0 = 1000.0
  
  sigmantr_out = zeros(SA.shape)
  
  # cycle through the depth index
  for k in range(len(z)):
  
    pt = pt_from_CT(CT[k, :, :], SA[k, :, :])   # convert CT to pt
    p  = p_from_z(z[k], lat)                    # convert z  to p
    ps = sp_from_SA(SA[k, :, :], p, lon, lat)   # convert SA to ps
    
    zws = sqrt( abs(ps) )
    
    # Numerator
    # T-Polynome
    zr1 = ( ( (-4.3159255086706703e-4 * pt + 8.1157118782170051e-2) * pt 
             + 2.2280832068441331e-1 ) * pt + 1002.3063688892480e0
           )
    # S-T Polynome
    zr2 = ( (-1.7052298331414675e-7 * ps - 3.1710675488863952e-3 * pt 
          - 1.0304537539692924e-4) * ps
          )
    # Denominator
    # T-Polynome
    zr3 = ( ( ( (-2.3850178558212048e-9 * pt -1.6212552470310961e-7) * pt
          + 7.8717799560577725e-5 ) * pt + 4.3907692647825900e-5  ) * pt + 1.0e0
          )
    # S-T Polynome
    zr4 = ( ( ( -2.2744455733317707e-9 * pt * pt + 6.0399864718597388e-6) * pt 
          - 5.1268124398160734e-4 ) * ps
          )
    # S-T Polynome
    zr5 = ( -1.3409379420216683e-9 * pt * pt - 3.6138532339703262e-5) * ps * zws

    # Neutral density
    sigmantr_out[k, :, :] = ( zr1 + zr2 ) / ( zr3 + zr4 + zr5 ) - zrau0
  
  return sigmantr_out

def pt_from_CT(CT, SA):
  """
  Compute POTENTIAL temperature from Conservative Temperature and 
  Absolute Salinity
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    CT           = Conservative Temperature      t        deg celsius
    SA           = Absolute Salinity             s        g / kg
    
  Returns:
    pt_out     = potential temperature         t        deg celsius
    
  """

  # ensures that SA is non-negative.
  SA = abs(SA)

  s1 = SA*0.995306702338459   # Note that 0.995306702338459 = (35./35.16504) 

  a0 = -1.446013646344788e-2
  a1 = -3.305308995852924e-3
  a2 =  1.062415929128982e-4
  a3 =  9.477566673794488e-1
  a4 =  2.166591947736613e-3
  a5 =  3.828842955039902e-3

  b0 =  1
  b1 =  6.506097115635800e-4
  b2 =  3.830289486850898e-3
  b3 =  1.247811760368034e-6

  a5CT = a5 * CT
  b3CT = b3 * CT
  CT_factor   = (a3 + a4 * s1 + a5CT)
  pt_num    = a0 + s1 * (a1 + a2 * s1) + CT * CT_factor
  pt_recden = 1.0 / ( b0 + b1 * s1 + CT * (b2 + b3CT) )
  pt        = pt_num * pt_recden      
                               # At this point the abs max error is 1.5e-2 deg C

  dpt_dCT = ( CT_factor + a5CT - (b2 + b3CT + b3CT) * pt ) * pt_recden
                            
  # start the 1.5 iterations through the modified Newton-Rapshon iterative 
  # method (McDougall and Wotherspoon, 2014). 

  CT_diff = CT_from_pt(pt, SA) - CT
  pt_old = pt
  pt = pt_old - CT_diff * dpt_dCT    # 1/2-way through the 1st modified N-R loop
                               # At this point the abs max error is 6.6e-5 deg C

  ptm = 0.5 * (pt + pt_old);

  # This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative 
  # of the Gibbs function with respect to temperature at zero sea pressure.  
  
  dpt_dCT = -3991.86795711963 / ( (ptm + 273.15) * gibbs_pt0_pt0(ptm, SA) )
  pt = pt_old - CT_diff * dpt_dCT       # end of 1st full modified N-R iteration
                               # At this point the abs max error is 1.e-10 deg C

  CT_diff = CT_from_pt(pt, SA) - CT
  pt_old = pt
  pt_out = pt_old - CT_diff * dpt_dCT
                                     # 1.5 iterations of the modified N-R method

  # The abs max error of the result is 1.42e-14 deg C
  
  return pt_out

################################################################################
# CT_from_pt (was gsw_CT_from_pt)                   CT --> potential temperature
#==========================================================================
#
# USAGE:
#  CT = gsw_CT_from_pt(SA,pt)
#
# DESCRIPTION:
#  Calculates Conservative Temperature of seawater from potential 
#  temperature (whose reference sea pressure is zero dbar).
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  pt  =  potential temperature (ITS-90)                          [ deg C ]
#
#  SA & pt need to have the same dimensions.
#
# OUTPUT:
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#
# AUTHOR: 
#  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
#  
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#    See section 3.3 of this TEOS-10 Manual. 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def CT_from_pt(pt, SA):
  """
  Compute CONSERVATIVE temperature from potential emperature and 
  Absolute Salinity
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    pt           = potential temperature         t        deg celsius
    SA           = Absolute Salinity             s          g / kg
    
  Returns:
    CT_out       = Conservative Temperature      t        deg celsius
    
  """

  # ensures that SA is non-negative.
  SA = abs(SA)

  sfac = 0.0248826675584615;                  # sfac = 1/(40.*(35.16504/35))

  x2 = sfac * SA
  x = sqrt(x2)
  y = pt * 0.025                           # normalize for F03 and F08.
  
  #--------------------------------------------------------------------------
  # The below polynomial for pot_enthalpy is the full expression for 
  # potential entahlpy in terms of SA and pt, obtained from the Gibbs 
  # function as below.  The above polynomial has simply collected like powers
  # of x and y so that it is computationally faster than calling the Gibbs 
  # function twice as is done in the commented code below.  When this code 
  # below is run, the results are identical to calculating pot_enthalpy as 
  # above, to machine precision.  
  #
  #  pr0 = zeros(size(SA));
  #  pot_enthalpy = gsw_gibbs(0,0,0,SA,pt,pr0) - ...
  #                       (273.15 + pt).*gsw_gibbs(0,1,0,SA,pt,pr0);
  #
  #-----------------This is the end of the alternative code------------------

  pot_enthalpy =  61.01362420681071 + (
      y * (168776.46138048015 +
      y * (-2735.2785605119625 + y * (2574.2164453821433 +
       y * (-1536.6644434977543 + y * (545.7340497931629 +
      (-50.91091728474331 - 18.30489878927802 * y) * y))))) +
      x2 * (268.5520265845071 + y * (-12019.028203559312 +
      y * (3734.858026725145 + y * (-2046.7671145057618 +
      y * (465.28655623826234 + (-0.6370820302376359 -
      10.650848542359153 * y) * y)))) +
      x * (937.2099110620707 + y * (588.1802812170108 +
      y * (248.39476522971285 + (-3.871557904936333 -
      2.6268019854268356 * y) * y)) + x * (-1687.914374187449 + 
      x * (246.9598888781377 +
      x * (123.59576582457964 - 48.5891069025409 * x)) +
      y * (936.3206544460336 +
      y * (-942.7827304544439 + y * (369.4389437509002 +
      (-33.83664947895248 - 9.987880382780322 * y) * y))))))
                                       )
      
  CT_out = pot_enthalpy / 3991.86795711963 # gsw_cp0 = 3991.86795711963

  return CT_out
