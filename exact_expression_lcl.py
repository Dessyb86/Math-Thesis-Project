def lcl(p,T,rh=None,rhl=None,rhs=None,return_ldl=False,return_min_lcl_ldl=False):

   import math
   import scipy.special

   # Parameters
   Ttrip = 273.16     # K
   ptrip = 611.65     # Pa
   E0v   = 2.3740e6   # J/kg
   E0s   = 0.3337e6   # J/kg
   ggr   = 9.81       # m/s^2
   rgasa = 287.04     # J/kg/K 
   rgasv = 461        # J/kg/K 
   cva   = 719        # J/kg/K
   cvv   = 1418       # J/kg/K 
   cvl   = 4119       # J/kg/K 
   cvs   = 1861       # J/kg/K 
   cpa   = cva + rgasa
   cpv   = cvv + rgasv

   # The saturation vapor pressure over liquid water
   def pvstarl(T):
      return ptrip * (T/Ttrip)**((cpv-cvl)/rgasv) * \
         np.exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) )
   
   # The saturation vapor pressure over solid ice
   def pvstars(T):
      return ptrip * (T/Ttrip)**((cpv-cvs)/rgasv) * \
         np.exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) )

   # Calculate pv from rh, rhl, or rhs
   rh_counter = 0
   if rh  is not None:
      rh_counter = rh_counter + 1
   if rhl is not None:
      rh_counter = rh_counter + 1
   if rhs is not None:
      rh_counter = rh_counter + 1
   if rh_counter != 1:
      print(rh_counter)
      exit('Error in lcl: Exactly one of rh, rhl, and rhs must be specified')
   if rh is not None:
      # The variable rh is assumed to be 
      # with respect to liquid if T > Ttrip and 
      # with respect to solid if T < Ttrip
      if T > Ttrip:
         pv = rh * pvstarl(T)
      else:
         pv = rh * pvstars(T)
      rhl = pv / pvstarl(T)
      rhs = pv / pvstars(T)
   elif rhl is not None:
      pv = rhl * pvstarl(T)
      rhs = pv / pvstars(T)
      if T > Ttrip:
         rh = rhl
      else:
         rh = rhs
   elif rhs is not None:
      pv = rhs * pvstars(T)
      rhl = pv / pvstarl(T)
      if T > Ttrip:
         rh = rhl
      else:
         rh = rhs
   if pv > p:
      return float('nan')

   # Calculate lcl_liquid and lcl_solid
   qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
   rgasm = (1-qv)*rgasa + qv*rgasv
   cpm = (1-qv)*cpa + qv*cpv
   if rh == 0:
      return cpm*T/ggr
   aL = -(cpv-cvl)/rgasv + cpm/rgasm
   bL = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
   cL = pv/pvstarl(T)*np.exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
   aS = -(cpv-cvs)/rgasv + cpm/rgasm
   bS = -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
   cS = pv/pvstars(T)*np.exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
   lcl = cpm*T/ggr*( 1 - \
      bL/(aL*scipy.special.lambertw(bL/aL*cL**(1/aL),-1).real) )
   ldl = cpm*T/ggr*( 1 - \
      bS/(aS*scipy.special.lambertw(bS/aS*cS**(1/aS),-1).real) )

   # Return either lcl or ldl
   if return_ldl and return_min_lcl_ldl:
      exit('return_ldl and return_min_lcl_ldl cannot both be true')
   elif return_ldl:
      return ldl
   elif return_min_lcl_ldl:
      return min(lcl,ldl)
   else:
      return lcl
