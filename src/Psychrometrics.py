'''

Pychrometrics - Imperial (English Units)

Library of typical Psychrometric functions 

** Dependencies **
1. Conversions.py - library of functions to convert between units

init development - 2017-09-30
updates - 2024-01-07

@author: Gregory Evans
@version: 1.0.0

'''


from PyChrometrics.src.Conversions import UnitConversions

from PyChrometrics.src.ThermodynamicProperties import MoistAirProperties

import pint

from math import exp
from math import log
from math import atan

class Pychrometrics():
    '''
    '''
    qty = pint.UnitRegistry.Quantity        # qty points to pint Quantity comprised of magnitude and unit
    u = pint.UnitRegistry()
    
    __properties = MoistAirProperties()
    convert = UnitConversions()
    __STD_ATM_PRESSURE = qty(14.696, 'psi')  # psia - std atm pressure at sea level (0 ft)
    __MIN_ELEV = qty(-1_400, 'ft')           # lowest acceptable elevation - Dead Sea (ft)  
    __MAX_ELEV = qty(30_000, 'ft')           # highest acceptable elevation - Mt. Everest (ft)
    
    
    def __init__(self,
                 input_val : qty,
                 input_type : str = "elevation"):
        '''
        @TODO - throw an Exception when either elevation or pressure is not provided
        '''
        self.__elevation = None
        self.__abs_pressure = None

        if input_type.lower() == "elevation":
            self.__set_elevation(input_val)
        elif input_type.lower() == "pressure":
            self.__set_pressure(input_val)


    def __set_elevation(self,
                        elevation : qty):
        
        if self.__MIN_ELEV < elevation <= self.__MAX_ELEV:
            self.__elevation = elevation
            self.__abs_pressure = self.P_atm_std(elevation)
        else:
            raise ValueError(f"Elevation is out of bounds. {self.__MIN_ELEV} < elevation <= {self.__MAX_ELEV}")
        
    def __set_pressure(self,
                       pressure: qty):
        
        MIN_PRES = self.P_atm_std(self.__MIN_ELEV)
        MAX_PRES = self.P_atm_std(self.__MAX_ELEV)

        if MAX_PRES <= pressure <= MAX_PRES:
            self.__abs_pressure = pressure
        else:
            raise ValueError(f"Pressure is out of bounds. {MIN_PRES} < pressure <= {MAX_PRES}")
        

    def P_atm_std(self,
                  elevation : qty) -> float:
        '''
        return: atmospheric pressure at elevation in psia
        '''
        k = self.qty(6.8754e-06, 'ft**-1')   # define const for calc of elevation

        P_atm = self.__STD_ATM_PRESSURE * (1 - elevation * k) ** 5.2559

        return P_atm
    
    def P_v_partial(self,
                    T_db : qty,
                    RH : qty ) -> qty:
        '''
        P_v_partial - vapor partial pressure
        '''
        P_ws = self.P_v_sat( T_db )
        
        return( RH * P_ws )
        
    def P_v_sat(self, 
                T_db : qty) -> qty:
        '''
        P_v_sat - saturated vapor partial pressure
        '''
               
        return( self.P_ice_sat( T_db ) if T_db < self.__properties.FREEZING_POINT() else self.P_w_sat( T_db ) )
    
    def P_ice_sat( self, 
                   T_db : qty ) -> qty:
        '''
        P_ice_sat - returns the water vapor saturation pressure over ice (psia)
        Hyland-Wexler Correlations - 1983 - ASHRAE 2001
        
        @param T_db - dry bulb temperature
        @return saturation vapor pressure over ice (psia)
        '''
                
        # m = [ -0.56745359e04,
        #       -.51523058,
        #       -0.009677843,
        #        0.62215701e-6,
        #        0.20747825e-08,
        #       -0.94840240e-12,
        #        0.41635019e01]

        C = [-5.6745359e03,
             6.3925247,
             -9.677843e-03,
             6.2215701e-07,
             2.0747825e-09,
             -9.4840240e-13,
             4.1635019e0 ]
        
        # convert dry bulb temperature in deg F to thermodynamic temperature in K
        T_abs = T_db.to('K')
        T = T_abs.magnitude

        ln_P_is = C[0] / T + C[1] + C[2] * T + C[3] * T**2 + C[4] * T**3 + C[5] * T**4 + C[6] * log(T)

        P_is = exp(ln_P_is)

        P_is = self.qty(P_is, 'Pa')

        return( P_is.to('psi') )
        
    def P_w_sat( self, 
                 T_db : qty) -> qty:
        '''
        P_w_sat - returns the water vapor saturation pressure over water (psia)
        Hyland-Wexler Correlations - 1983 - ASHRAE 2001
        
        @param T_db - Thermodynamic temperature ( K )
        @return saturation vapor pressure over water (psia)
        '''
        
        # h = [ -0.58002206e04,
        #       -5.516256,
        #       -0.48640239e-01,
        #        0.41764768e-04,
        #       -0.14452093e-07,
        #        6.5459673e00 ]
        
        C = [ -5.8002206e03,
              1.3914994e00,
              -4.8640239e-02,
               4.1764768e-05,
              -1.4452093e-08,
               6.5459673e00 ]
        
        # convert dry bulb temperature in deg F to thermodynamic temperature in K
        T_abs = T_db.to('K')
        T = T_abs.magnitude

        ln_P_sat = C[0] / T + C[1] + C[2] * T + C[3] * T**2 + C[4] * T**3 + C[5] * log(T)

        P_sat = exp(ln_P_sat)

        P_sat = self.qty(P_sat, 'Pa')

        return( P_sat.to('psi') )
    
    def T_wb_iter(self, 
                  T_db : float,
                  RH : float) -> float:
        ''' 
        wet bulb temperature - iterative method 
        
        Adapted from research listed in:
        Al-Ismaili, A.M., Al-Azri, N.A. (2016). Simple Iterative Approach to Calculate
        Wet-Bulb Temperature for Estimating Evaporative Cooling Efficiency
        International Journal of Agriculture Innovations and Research, 4(6), (Online) 1013-1018.
        http://ijair.org/administrator/components/com_jresearch/files/publications/IJAIR_1945_Final.pdf
        
        @param - T_db: dry bulb Temperature (deg F)
        @param - RH: relative humidity (fractional form)
        @return - wet bulb Temperature (deg F)
        
        @TODO - complete implementation of method with W0 / W1 methods
        '''

        T_wb = T_db                     # Initial guess at Twb
        W = self.W( T_db, RH )
        W0 = self.W0( T_db, T_wb )
        epsilon = 1e-6                  # stop condition
        
        def deltaW():
            return (( abs( W - W0) / W )/ 100 )
        
        while ( abs( W0 - W ) > epsilon ):
            
            T_wb *= 1 - deltaW() if abs( W0 - W ) > 0 else 1 + deltaW()
            
            W0 = self.W0( T_db, T_wb )
        
        return T_wb
    
    def T_wb_regression(self, 
                        T_db : float,
                        RH : float) -> float:
        ''' 
        wet bulb temperature - regression method
        
        Adapted from research listed in:
        Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature. 
        Journal of Applied Meteorology and Climatology, 50(11), (Online) 2267-2269. doi:10.1175/jamc-d-11-0143.1 
        http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
        
        @param - dry bulb Temperature (deg F)
        @param - relative humidity (fractional form)
        @return - wet bulb Temperature (deg F)
        
        '''
        
        T_db_c = self.convert.F_to_C( T_db )
        
        RH = self.convert.convert_RH( RH )
                
        # 'coefficients' of regression equation
        # Twb = Tdb * a( RH )  + b( Tdb, RH) - c( RH ) + d( RH ) - e
        
        a = atan( ( 0.151977 * ( RH + 8.313659 ) ** 0.5 ) )
        bT = atan( T_db + RH ) 
        c = atan( RH - 1.676331 )
        d = 0.00391838 * RH ** 1.5 * atan( 0.023101 * RH )
        e = 4.686035
        
        T_wb_c = T_db_c * a + bT - c + d - e
         
        # convert from Celsius to Farenheit
        T_wb = self.convert.C_to_F(( T_wb_c ) )
        
        return T_wb
    
    def T_dp( self, T_db, RH ):

        c = [ 100.45,
              33.193,
              2.319,
              0.17074,
              1.2063 ]
        
        P_w = self.P_v_partial( T_db, RH )       # psia
        a = log( P_w )                           # ln( psia )
        T_dp = 0                                 # deg C
        
        if ( T_db >= 32.0 ):
            for i in range( len( c ) ):
                T_dp += c[ i ] * a ** i if i < 4 else c[ i ] * P_w ** 0.1984
        else:
            T_dp += 90.12 + 26.142 * a + 0.8927 * a ** 2.0

        return ( T_dp )
        
    
    def RH( self, T_db, T_dp ):
        '''
        RH - relative humidity
        '''
        Tdb_cel = self.convert.F_to_C(T_db)
        Tdp_cel = self.convert.F_to_C(T_dp)

        a = 17.625
        b = 243.04

        A = exp( (a * Tdp_cel)/(b + Tdp_cel) )
        B = exp( (a * Tdb_cel)/(b + Tdb_cel) )

        return 100 * A/B
    
    def W(self,
          T_db : qty,
          RH : qty ):
        '''
        Humidity Ratio: W - function calculates the humidity ratio given dry bulb temperature and relative humidity
        
        @param - T_db: dry bulb temperature (deg F)
        @param - RH: relative humidity ( lbm vapor / lbm dry air )
        @return - humidity ratio (lbm/lbm)
        
        '''
        
        P_atm = self.__abs_pressure
        
        # Pws = self.P_v_sat( T_db )
        Pws = self.P_w_sat( T_db )
        
        MMR = self.__properties.MMR()
        
        return ( MMR * ( Pws * RH /( P_atm - Pws * RH) ) ) # equation is from ASHRAE 2002 / 2005
    
    def W_sat( self, T_db):
        '''
        saturated humidity ratio: sat_W - function calculates the saturated humidity ratio based on dry bulb temperature
        
        @param - T_db: dry bulb temperature
        @return - saturated humidity ratio (lbm/lbm)
        
        '''
        
        P_atm = self.__abs_pressure
        
        Pws = self.P_v_sat( T_db )
        
        MMR = self.__properties.MMR()
        
        return ( MMR * Pws /( P_atm - Pws ) ) # equation is from ASHRAE 2002 / 2005
    
    def W0 ( self, Tdb, Twb ):
        
        '''
        W0 - humidity ratio calculation as a function of Tdb and Twb with variable specific heats for air, water, vapor
        
        Adaped from research listed in:
        Al-Ismaili, A.M., Al-Azri, N.A. (2016). Simple Iterative Approach to Calculate
        Wet-Bulb Temperature for Estimating Evaporative Cooling Efficiency
        International Journal of Agriculture Innovations and Research, 4(6), (Online) 1013-1018.
        http://ijair.org/administrator/components/com_jresearch/files/publications/IJAIR_1945_Final.pdf
        
        @param - T_db: dry bulb temperature in degrees F
        @param - T_wb: wet bulb temperature in degrees F
        @return - humidity ratio based on constant specific heat of air, water, and vapor (lbm/lbm)
        
        '''
    
        h_fg = self.__properties.h_fg()                # Latent Heat of Vaporization (Btu/lbm)
        
        Cp_w = self.__properties.Cp_water                                   # Cp of water function
        Cp_da = self.__properties.Cp_dry_air                                # Cp of dry air function
        Cp_v = self.__properties.Cp_vapor                                   # Cp of water vapor function
        Ws_wb = self.W_sat                                     # saturated humidity ratio Ws(Twb)  
        
        a = ( ( h_fg - ( Cp_w( Twb ) - Cp_v( Tdb ) ) * Twb ) * Ws_wb( Twb ) - Cp_da( Tdb, Twb ) * ( Tdb - Twb ) )
        b = ( h_fg + Cp_v( Tdb ) * Tdb - Cp_w( Twb ) * Twb )
        
        return ( a / b )
    
    def W1 ( self, T_db, T_wb ):
        
        '''
        humidity ratio - calculation as a function of Tdb and Twb with const specific heats for air, water, vapor
        
        Adaped from research listed in:
        Al-Ismaili, A.M., Al-Azri, N.A. (2016). Simple Iterative Approach to Calculate
        Wet-Bulb Temperature for Estimating Evaporative Cooling Efficiency
        International Journal of Agriculture Innovations and Research, 4(6), (Online) 1013-1018.
        http://ijair.org/administrator/components/com_jresearch/files/publications/IJAIR_1945_Final.pdf
        
        @param - T_db: dry bulb temperature in degrees F
        @param - T_wb: wet bulb temperature in degrees F
        @return - humidity ratio based on constant specific heat of air, water, and vapor (lbm/lbm)
        
        '''
        
        h_fg = self.__properties.h_fg()                    # Latent Heat of Vaporization (Btu/lbm)        
        Cp_da = self.__properties.Cp_da_const()                        # specific heat of dry air
        Cp_w = self.__properties.Cp_w_const()                           #  specific heat of water 
        Cp_v = self.__properties.Cp_v_const()                           # specific heat of vapor 
        Ws_wb = self.W_sat                                         # saturated humidity ratio function (pointer) Ws(Twb)  
        
        a = ( ( h_fg - ( Cp_w - Cp_v ) * T_wb ) * Ws_wb( T_wb ) - Cp_da * ( T_db - T_wb ) )
        b = ( h_fg + Cp_v * T_db - Cp_w * T_wb )
        
        return ( a / b )
    
    def grains_moisture( self, T_db, RH ):
        ''' 
        Returns the grains of moisture 
        '''
    
        W = self.W( T_db, RH )
        grains = self.__properties.GRAINS_PER_LBM()

        return ( grains * W )
    
    def specific_humidity( self ):
        pass
    
    def absolute_humidity( self ):
        pass
    
    def degree_saturation( self, T_db, RH ):
        '''
        Returns the degree of saturation
        '''
        
        W = self.W( T_db, RH )
        Ws = self.W_sat( T_db )
        
        return ( W / Ws * 100)
    
    def specific_volume( self, T_db, RH ):
        '''
        '''
        T_db_R = self.convert.F_to_R( T_db )
        W = self.W(T_db, RH)
        P = self.__abs_pressure
        
        return 0.3704 * T_db_R * ( 1 + 1.6078 * W ) / P
    
    def density( self, T_db, RH ):
        '''
        Returns the density of moist air based on the specific volume
        '''
        v = self.specific_volume(T_db, RH)
        rho = 1 / v
        return rho
    
    def enthalpy_dry_air( self, T_db, RH ):
        '''
        Returns the enthalpy of dry air
        h_da = Cp_da * T_db
        '''
    
        T_wb = self.T_wb_iter(T_db, RH)
        Cp_da = self.__properties.Cp_dry_air(T_db, T_wb)
        
        return( Cp_da * T_db )
    
    def enthalpy_vapor( self, T_db, RH ):
        '''
        Returns the latent enthalpy
        h_v = W * (c_pv * T + h_fg)
        '''

        W = self.W( T_db, RH )
        C_pv = self.__properties.Cp_vapor( T_db )
        h_fg = self.__properties.h_fg()

        return W * (C_pv * T_db + h_fg)
    
    def enthalpy( self, T_db, RH ):
        '''
        Returns the enthalpy of moist air (latent and sensible)
        h = h_da + h_v
        '''
        h_da = self.enthalpy_dry_air(T_db, RH)
        h_v = self.enthalpy_vapor(T_db, RH)

        return h_da + h_v

    def pressure(self):
        '''
        Returns the absolute pressure
        '''
        return( self.__abs_pressure )





