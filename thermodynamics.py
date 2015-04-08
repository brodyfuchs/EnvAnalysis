from numpy import *

"""
A series of functions to calculate thermodynamic variables

compiled by Brett Basarab
Colorado State University
bbasarab[at] atmos.colostate.edu
last updated March 2015

Unless otherwise indicated, most of these functions were written
by Thomas Chubb and downloaded as part of a free python skew-T plotting
package

Most functions take arguments in standard SI units:
    temperature in K
    pressure in Pa
    mixing ratio in kg/kg
Some take temperature in Celsius
     
"""

#-----------------------------------------------------------------------
# Here we go. A set of functions that I use from time to time to calculate 
# the basic stuff that I'm sick of doing over and over! I'm going to 
# endeavour to include references and global constants to make it all nice 
# and legible.
#-----------------------------------------------------------------------

Rs_da=287.05          # Specific gas const for dry air, J/kg/K
Rs_v=461.51           # Specific gas const for water vapour, J/kg/K
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218	      # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)

_devel="working"


def Cel2K(tempc):
    return tempc+273.15


def Cel2F(tempc):
    return tempc*(9/5.)+32.


def F2K(tempf):
    return 273.15+(5/9.)*(tempf-32.) 


def F2Cel(tempf):
    return (5/9.)*(tempf-32.)


def Theta(tempk,pres,pref=100000.):
    """Potential Temperature

    INPUTS: 
    tempk (K)
    pres (Pa)
    pref: Reference pressure (default 100000 Pa)

    OUTPUTS: Theta (K)

    Source: Wikipedia
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """

    try:
	minpres=min(pres)
    except TypeError:
	minpres=pres

    if minpres<2000:
	print "WARNING: P<2000 Pa; did you input a value in hPa?"

    return tempk*(pref/pres)**(Rs_da/Cp_da)

def TempK(theta,pres,pref=100000.):
    """Inverts Theta function. i.e, gets temperature at pres
    based on theta."""

    try:
	minpres=min(pres)
    except TypeError:
	minpres=pres

    if minpres<2000:
	print "WARNING: P<2000 Pa; did you input a value in hPa?"

    return theta*(pres/pref)**(Rs_da/Cp_da)

def ThetaE_R():
    """Equivalent potential temperature:
    Exact equation for reversible processes"""
    raise NotImplementedError
  
 
def ThetaE_Pseudo():
    """Equivalent potential temperature:
    Empirical formula valid for pseudo adiabatic processes
    from Bolton (1980)"""
    raise NotImplementedError

 
def ThetaE_App(tempk,pres,w,pref=100000.):
    """Equivalent potential temperature: approximate formula"""
    theta = Theta(tempk,pres,pref=pref)
    theta_e = theta*exp(Lv*w/(Cp_da*tempk))

    return theta_e


def ThetaES(tempk,pres,pref=100000.):
    """Saturated equivalent potential temperature: approximate formula"""
    theta = Theta(tempk,pres,pref=pref)
    ws = SatMixRatio(tempk,pres)
    theta_es = theta*exp(Lv*ws/(Cp_da*tempk))

    return theta_es


def ThetaV(tempk,pres,e):
    """Virtual Potential Temperature
    
    INPUTS
    tempk (K)
    pres (Pa)
    e: Water vapour pressure (Pa) (Optional)
    """ 

    mixr=MixRatio(e,pres)
    theta=Theta(tempk,pres)

    return theta*(1+mixr/Epsilon)/(1+mixr)


def DSE(tempk,height):
    """Dry static energy"""
    return Cp_da*tempk+grav*height


def MSE(tempk,height,w):
    """Moist static energy"""
    return Cp_da*tempk+grav*height+Lv*w


def SMSE(tempk,height,pres):
    """Saturation moist static energy"""
    ws = SatMixRatio(tempk,pres) # approximate formula
    
    return Cp_da*tempk+grav*height+Lv*ws


def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """
    
    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=MixRatio(es,pres)

    if e is None:
	# assume saturated
	e=es

    w=MixRatio(e,pres)
    #tempk = tempc+degCtoK
    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma


def DensMoist(tempk,pres,mixr):
    """Density of moist air"""
    
    virtualT=VirtualTempFromMixR(tempk,mixr)
    return pres/(Rs_da*virtualT)


def VirtualTemp(tempk,pres,e):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    e: vapour pressure (Pa)
    p: static pressure (Pa)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia)."""

    tempvk=tempk/(1-(e/pres)*(1-Epsilon))
    return tempvk
    

def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk*(1.0+0.6*mixr)


def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """

    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)


def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice

   
    RETURNS: e_sat  (Pa)
    
    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)
    
    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """
    
    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
	# return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
	return over_liquid
    elif phase=="ice":
	# return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
	return where(tempc<0,over_ice,over_liquid)
    else:
	raise NotImplementedError


def MixRatio(e,pres):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    pres (Pa) Ambient pressure
          
    RETURNS
    w (kg/kg) Water vapor mixing ratio`
    """

    return Epsilon*e/(pres-e)


def SatMixRatio(tempk,pres):
    """Calculate saturation mixing ratio of water
    vapor with respect to liquid water, given 
    temperature and pressure"""

    tempc = tempk-degCtoK
    es = SatVap(tempc)
 
    return Epsilon*es/(pres-es) # approximate formula


def MixR2VaporPress(w,pres):
    """Return Vapor Pressure given Mixing Ratio and Pressure
    INPUTS
    w (kg kg^-1) Water vapor mixing ratio`
    pres (Pa) Ambient pressure
          
    RETURNS
    e (Pa) Water vapor pressure
    """

    return w*pres/(Epsilon+w)


def MixR2Q(w):
    """Simple conversion from mixing ratio to specific
    humidity so I don't forget the difference between them"""
    q = w/(w+1.)
    
    return q


def Dwpt2VapPres(dwpt):
    """Water vapor pressure
    INPUTS
    dwpt (C) Dew Point Temperature (for SATURATION vapor 
	     pressure use tempk)
          
    RETURNS
    e (Pa) Water Vapor Pressure

    SOURCE:
    Bolton, Monthly Weather Review, 1980, p 1047, eq. (10)
    """

    return 611.2*exp(17.67*dwpt/(243.5+dwpt))


def VapPres2Dwpt(e):
    """ Use Bolton's (1980, MWR, p1047) formulae to find tdew.
    INPUTS:
    e (Pa) Water Vapor Pressure
    OUTPUTS:
    Td (C) 
      """

    ln_ratio=log(e/611.2)
    Td=((17.67-ln_ratio)*degCtoK+243.5*ln_ratio)/(17.67-ln_ratio)
    return Td-degCtoK

def interp_parcel_path(pres_s,pres,temp,height,p,tdry,pwet,twet):
    """
    Interpolates the output from lift_parcel in the SkewT module in order to calculate
    a parcel temperature at each vertical level of your domain. <pres>, <temp>, and <height>
    may go upward from surface (so <pres>,<temp> decrease, <height> increases) or downward from
    top of domain. (But make sure all three of these arrays go in the same direction). If pressure
    array was input in Pa units, will convert to hPa. 

    Arguments (all are required):
	pres_s: surface pressure (hPa) 
	pres: your pressure array (hPa)
	temp: your temperature array (Celsius)
        height: your height array (meters)
        p: pressures along dry adiabat output by SkewT.lift_parcel (hPa) 
        tdry: temperatures along dry adiabat output by SkewT.lift_parcel (Celsius)
        pwet: pressures along moist adiabat output by SkewT.lift_parcel (hPa)
        twet: temperatures along moist adiabat output by SkewT.lift_parcel (Celsius)

    Returns: 
        parcel_temp: array of parcel temperatures interpolated to your vertical grid;
        will go from top of domain to surface
        h_el,p_el: height (meters) and pressure (hPa) of the equilibrium level (EL)
        h_lfc,p_lfc: height (meters) and pressure (hPa) of the level of free convection (LFC)
        h_frz,p_frz: height (meters) and pressure (hPa) of the environmental freezing level
        h_frz_parc,p_frz_parc: height (meters) and pressure (hPa) of the parcel path freezing level
        e,l,frz: indices of your arrays corresponding to the EL, LFC, and environmental freezing level
    """
    
    # the calculations below operate on arrays going from top of domain to surface
    # so check the direction of input arrays (<pres>,<temp>,<height>) here
    if all(diff(pres)<=0): # arrays start at surface so flip them
        pres = pres[::-1]
        temp = temp[::-1]
        height = height[::-1]

    if pres.max()>1500.:
       print 'MAX PRESSURE IS > 1500. ASSUMING PRESSURE WAS INPUT IN PA, CONVERTING TO HPA'
       pres = pres/100.
    
    # because of how pwet and twet are set up, exlude all indices above 100-hPa
    cutoff_pres = pres>=100.
    pres = pres[cutoff_pres]
    temp = temp[cutoff_pres]
    height = height[cutoff_pres]
    p_lcl = pwet[0] # pwet[0] and twet[0] are pressure and temperature of the LCL
    h_lcl = interp(p_lcl,pres,height)
    print 'p_lcl, h_lcl:', p_lcl, h_lcl

    atm = pres <= pres_s # finding pres elements that are above the surface
    p_atm = pres[atm]
    pres = pres[atm]
    temp = temp[atm]
    height = height[atm]
    
    dry = where(p_atm >= p_lcl) # indices of pressure levels below LCL (higher pressure than LCL)
    moist = where(p_atm < p_lcl) # indices of pressure levels above LCL
 
    parcel_temp = zeros(p_atm.shape[0]) # (degrees C); levels below surface will be left at zero
 
    # below LCL (parcel follows dry adiabat)
    # interpolate output from lift_parcel to your model or observational levels
    for d in dry[0]: # assign dry temperatures from LCL downward
	p_diff = p-p_atm[d] # 100-element pressure array minus element of interest in your pres array
	below = where(p_diff >= 0)[0][-1] # index of closest pressure in p level just below the LCL
	above = below + 1 # index of closest pressure level in p just above LCL
	# parcel temp follows dry adiabat below LCL
        # get dry parcel temperatures by interpolating on the function tdry(p)
        # p and tdry are decreasing, np.interp only works on increasing arrays (hence the above then below ordering)
	parcel_temp[d] = interp(p_atm[d], [p[above], p[below]], [tdry[above], tdry[below]])

    # can't interpolate on the end; parcel_temp is an increasing array (goes from top of domain), 
    # so just make the first element be the highest-up wet parcel path temperature
    # OK that parcel_temp is increasing because so is temp 
    parcel_temp[0] = twet[-1] # parcel_temp increases, and can't interpolate on the end;

    # above LCL (parcel follows moist adiabat)
    for m in moist[0][1:]: # assign wet temperatures from top of domain downward
	p_diff = pwet-p_atm[m]
	below = where(p_diff >= 0)[0][-1] # closest pressure level just below the LCL
	above = below + 1
	# parcel temp follows moist adiabat above LCL
	parcel_temp[m] = interp(p_atm[m], [pwet[above], pwet[below]], [twet[above], twet[below]])

    p_el,h_el,p_lfc,h_lfc = 0,0,0,0 # initialize these variables
    # find the equilibrium level (EL)
    for i in range(temp.shape[0]):
	t_diff = temp[i] - parcel_temp[i]
        #print temp[i],parcel_temp[i],'t_diff: ',t_diff
	# this condition seems backwards for finding the EL, but remember that temp and parcel_temp start at
	# TOP of domain. So you're finding FIRST instance of parcel being warmer or equal temperature to environment,
	if t_diff <= 0:
            print 'FOUND EL'
	    below = i
	    e = i # index of EL
	    break # break out of the loop
        else: # haven't found EL, means EL, LFC, and therefore CAPE quantities are undefined
            below = i
            e = i
         
    if e==temp.shape[0]-1: # didn't find an EL
	l = e # set LFC index to EL index so that CAPE calculation will be zero
	p_el,h_el,p_lfc,h_lfc = -1,-1,-1,-1 # set these values to undefined
	print 'EL NOT FOUND'

    #print 'for EL: below, e: ',below,e

    if (p_el!=-1) and (h_el!=-1): # MEANS THAT EL HAS BEEN FOUND
	# interpolate on function P(T-T_parcel); find pressure where temp-parcel_temp = 0
	# if t_diff = 0, then the interpolation is trivial, below is the index of your EL
	# if t_diff < 0, then you just a little bit below the EL; interpolate between (temp[below]-parcel_temp[below]),
	# and (temp[below-1]-parcel_temp[below-1]) 
	p_el = interp(0, [temp[below]-parcel_temp[below], temp[below-1]-parcel_temp[below-1]],
	       [p_atm[below],p_atm[below-1]]) # pressure of equilibrium level
	#print 'heights near EL: ', height[below],height[below-1]
	h_el = interp(0, [temp[below]-parcel_temp[below], temp[below-1]-parcel_temp[below-1]],
	       [height[below],height[below-1]]) # height of equilibrium level

	if p_el<p_lcl: #### MEANS THAT EL IS PHYSICAL 
	    ### CALCULATE PRESSURE AND HEIGHT OF THE LEVEL OF FREE CONVECTION
	    for i in range(below, temp.shape[0]):
		t_diff = temp[i] - parcel_temp[i]
		#print 't_diff: ',t_diff
		# now, starting from the EL, move downward in height space to find the first instance of parcel temperature
		# being cooler or equal temperature to environment, i.e., the LFC
		found_lfc = 0
		if t_diff >= 0: # parcel is cooler than environment
		    #print 'found LFC'
		    found_lfc=1
		    below = i
		    l = i # index of LFC
		    break
		else: # haven't found the LFC, means parcel is warmer than env. all the way to surface
		    below = i
		    l = i
	    #print 'for LFC: below, l, p_atm.shape: ',below,l,p_atm.shape[0]
	    #print 'temp.shape, parcel_temp.shape: ',temp.shape[0],parcel_temp.shape[0]
	    #print 'parcel_temp: ',parcel_temp
	    if l == p_atm.shape[0]-1: # got to end of array, no LFC found 
		p_lfc = p_lcl
		h_lfc = h_lcl
	    #elif found_lfc == 0:
		#p_lfc = p_lcl
		#h_lfc = h_lcl
	    else: # 
		p_lfc = interp(0, [temp[below-1]-parcel_temp[below-1], temp[below]-parcel_temp[below]],
				    [p_atm[below-1],p_atm[below]])
	    #    h_lfc = interp(0, [temp[below-1]-parcel_temp[below-1], temp[below]-parcel_temp[below]], 
	    #                            [ht[below-1],ht[below]])
		h_lfc = interp(p_lfc, pres, height)
	else: #### EL IS NOT PHYSICAL
	    p_el,h_el,p_lfc,h_lfc = -1,-1,-1,-1 # set these values to undefined
	    l = e # set LFC index to EL index so that CAPE calculation will be zero

    print 'p_el, h_el: ',p_el,h_el
    print 'p_lfc, h_lfc: ',p_lfc,h_lfc

    # find height of the (environmental) freezing level
    frz = where(abs(temp) == abs(temp).min())[0][0]
    p_frz = pres[frz]
    h_frz = height[frz]
   
    # find the height of the (parcel path) freezing level
    i_frz_parc = where(abs(parcel_temp[pres<=pres_s]) == abs(parcel_temp[pres<=pres_s]).min())[0][0]
    h_frz_parc = height[i_frz_parc]
    p_frz_parc = pres[i_frz_parc]     

    #print 'DONE INTERPOLATING PARCEL PATH!!\n\n'
    # returns array of parcel temperatures, pressure and height of LFC,
    #pressure and height of EL, index of freezing height, LFC, and EL
    return parcel_temp,p_lcl,h_lcl,p_lfc,h_lfc,p_el,h_el,p_frz,h_frz,p_frz_parc,h_frz_parc,frz,l,e

def calc_cape(temp,height,parcel_temp,h_lfc,l,h_el,e,frz):
    """
    This function calculates the convective available potential energy (CAPE) between the level
    of free convection and the equilibrium level. It also calculates a "warm" CAPE value (CAPE
    below the freezing level). <temp>, <height>, and <parcel_temp> may go upward from surface 
    (so <temp>, <parcel_temp> decrease, <height> increases) or downward from top of domain. (But 
    make sure all three of these arrays go in the same direction). 

    Arguments: 
        temp: array of temperatures (Celsius)
        height: array of heights (meters)
        parcel_temp: array of parcel temperatures interpolated to <temp>, <height> grid (Celsius)
        h_lfc: height of the LFC (meters)
        l: index of <parcel_temp>-<temp> closest the LFC
        h_el: height of the EL (meters)
        e: index <parcel_temp>-<temp> closest to the EL
        frz: index closest to the freezing level
    """

    # the calculations below operate on arrays going from top of domain to surface
    # so check the direction of input arrays (<pres>,<temp>,<height>) here
    if all(diff(height)>0): # arrays start at surface so flip them
        temp = temp[::-1]
        height = height[::-1]
        parcel_temp = parcel_temp[::-1]

    ### CAPE CALCULATIONS ###
    cape = 0.
    warm_cape = 0.
    if l > e: # NEED TO MAKE SURE THERE IS A REASONABLE EL
        print 'l > e; calculating CAPE'
	for i in range(e,l-1): # descretized CAPE calculation: looping from EL to LFC to get CAPE 
	    cape += grav*(parcel_temp[i]-temp[i])*(height[i]-height[i+1])/(temp[i]+degCtoK)
	    if i >= frz: # CAPE below the freezing level woot!!
		warm_cape += grav*(parcel_temp[i]-temp[i])*(height[i]-height[i+1])/(temp[i]+degCtoK)
	# ADD ON A LITTLE NEAR THE EL
	el_cape = grav*(parcel_temp[e]-temp[e])*(h_el-height[e])/(temp[e]+degCtoK)
	# SUBTRACT OFF ANY BELOW THE LFC
	lfc_cape = grav*(parcel_temp[l]-temp[l])*(h_lfc-height[l])/(temp[l]+degCtoK)
        print 'EL CAPE, LFC CAPE:',el_cape,lfc_cape
	if (abs(el_cape) > 1000) | (abs(lfc_cape) > 1000):
	    print 'BAD near-EL or near-LFC CAPE values; not including'
	else:
	    cape += el_cape
	    cape += lfc_cape
	    warm_cape += lfc_cape
	del lfc_cape # delete this variable
	del el_cape
	if cape < 0: cape = 0
    
    if cape!=0.: 
	ncape = cape/(h_el-h_lfc)
    else:
        ncape = 0.    

    return cape,ncape,warm_cape
