'''

PsychroCalcs - Imperial (English Units)

Library of typical Psychrometric functions 

** Dependencies **
1. Conversions.py - library of functions to convert between units

init development - 30/09/2017

@author: Gregory Evans
@version: 1.0.0

'''


from Conversions import UnitConversions

from ThermodynamicProperties import MoistAirProperties

from math import exp
from math import log
from math import atan

class PsychroCalcs():
    '''
    '''
    
    __properties = MoistAirProperties()
    convert = UnitConversions()
    
    __elevation = 0                     # ft - default elevation of 0
    __abs_pressure = 0                  # psia - default absolute atmospheric pressure of 0 psia
    
    
    def __init__(
                 self,
                 elevation = None,
                 abs_pressure = None,
                 ):
        '''
        @TODO - calculate elevation given pressure
              - throw an Exception when either elevation or pressure is not provided
        '''
        
        self.__elevation = elevation if (elevation != None) else 0.0
        
        self.__abs_pressure = abs_pressure if abs_pressure != None else self.__std_atm_pressure()
           
    def __std_atm_pressure( self ):
        '''
        '''
        return ( 14.696 * ( 1 - self.__elevation * 6.8754e-06 ) ** 5.2559 )
    
    def __partial_vap_pressure( self, T_db, RH ):
        '''
        '''
        
        P_ws = self.__sat_vapor_pressure( T_db )
        
        return( RH * P_ws )
        
    def __sat_vapor_pressure( self, T_db ):
        
        return( self.__sat_pressure_ice( T_db ) if T_db < self.__properties.FREEZING_POINT() else self.__sat_pressure_water( T_db ) )
    
    def __sat_pressure_ice( self, T_db ):
        '''
        __sat_pressure_ice - returns the water vapor saturation pressure over ice (psia)
        Hyland-Wexler Correlations - 1983 - ASHRAE 2001
        
        @param T_db - dry bulb temperature
        @return saturation vapor pressure over ice (psia)
        '''
                
        m = [ -0.56745359e04,
              -.51523058,
              -0.009677843,
               0.62215701e-6,
               0.20747825e-08,
              -0.94840240e-12,
               0.41635019e01]
           
        # ln(Pw) = Sum(m_i * T^i-1, i = 0, 5) + m_6 * ln(T)
        # natural log of saturation vapor pressure for ice
        ln_P_is = 0
        
        # convert dry bulb temperature in deg F to thermodynamic temperature in K
        T_abs = self.convert.C_to_K( self.convert.F_to_C( T_db ) )
        
        for i in range( -1, len( m ) - 1 ):
            
            ln_P_is += m[ i + 1 ] * T_abs ** i if ( i < len( m ) - 2 ) else m[ i + 1 ] * log( T_abs )
                              
        P_is = exp( ln_P_is )
                    
        return ( self.convert.Pa_to_psia( P_is ) )
        
    def __sat_pressure_water( self, T_db ):
        '''
        __sat_pressure_water - returns the water vapor saturation pressure over water (psia)
        Hyland-Wexler Correlations - 1983 - ASHRAE 2001
        
        @param T_db - Thermodynamic temperature ( K )
        @return saturation vapor pressure over water (psia)
        '''
        
        h = [ -0.58002206e04,
              -5.516256,
              -0.48640239e-01,
               0.41764768e-04,
              -0.14452093e-07,
               6.5459673e00 ]

        # ln(Pws) = Sum(h_i * T^1, i = -1, 3) + h_r * ln(T)    
        # natural log of saturation vapor pressure for water
        ln_P_ws = 0
        
        # convert dry bulb temperature in deg F to thermodynamic temperature in K
        T_abs = self.convert.C_to_K( self.convert.F_to_C( T_db ) )
        
        for i in range( -1, len( h ) - 1 ):
            
            ln_P_ws += h[ i + 1 ] * T_abs ** i if ( i < len( h ) - 2 ) else h[ i + 1 ] * log( T_abs )
           
        P_ws = exp( ln_P_ws )

        return ( self.convert.Pa_to_psia( P_ws ) )
        
    def T_db( self ):
        pass
    
    def T_wb_iter( 
                  self, 
                  T_db,
                  RH 
                  ):
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
        W = self.humidity_ratio( T_db, RH )
        W0 = self.W0( T_db, T_wb )
        epsilon = 1e-6                  # stop condition
        
        def deltaW():
            return (( abs( W - W0) / W )/ 100 )
        
        while ( abs( W0 - W ) > epsilon ):
            
            T_wb *= 1 - deltaW() if abs( W0 - W ) > 0 else 1 + deltaW()
            
            W0 = self.W0( T_db, T_wb )
        
        return T_wb
    
    def T_wb_regression( 
                        self, 
                        T_db, 
                        RH 
                        ):
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
        
        P_w = self.__partial_vap_pressure( T_db, RH )       # psia
        a = log( P_w )                                      # ln( psia )
        T_dp = 0                                            # deg C
        
        if ( T_db >= 32.0 ):
            for i in range( len( c ) ):
                T_dp += c[ i ] * a ** i if i < 4 else c[ i ] * P_w ** 0.1984
        else:
            T_dp += 90.12 + 26.142 * a + 0.8927 * a ** 2.0

        return ( T_dp )
        
    
    def RH( self ):
        pass
    
    def humidity_ratio( self, T_db, RH ):
        '''
        humidity_ratio - function calculates the humidity ratio given dry bulb temperature and relative humidity
        
        @param - T_db: dry bulb temperature (deg F)
        @param - RH: relative humidity ( lbm vapor / lbm dry air )
        @return - humidity ratio (lbm/lbm)
        
        '''
        
        P_atm = self.__abs_pressure
        
        Pws = self.__sat_vapor_pressure( T_db )
        
        MMR = self.__properties.get_molec_mass_ratio()
        
        return ( MMR * ( Pws * RH /( P_atm - Pws * RH) ) ) # equation is from ASHRAE 2002 / 2005
    
    def sat_humidity_ratio( self, T_db):
        '''
        sat_humidity_ratio - function calculates the saturated humidity ratio based on dry bulb temperature
        
        @param - T_db: dry bulb temperature
        @return - saturated humidity ratio (lbm/lbm)
        
        '''
        
        P_atm = self.__abs_pressure
        
        Pws = self.__sat_vapor_pressure( T_db )
        
        MMR = self.__properties.get_molec_mass_ratio()
        
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
    
        h_fg = self.__properties.get_enthalpy_vaporization()                # Latent Heat of Vaporization (Btu/lbm)
        
        Cp_w = self.__properties.Cp_water                                   # Cp of water function
        Cp_da = self.__properties.Cp_dry_air                                # Cp of dry air function
        Cp_v = self.__properties.Cp_vapor                                   # Cp of water vapor function
        Ws_wb = self.sat_humidity_ratio                                     # saturated humidity ratio Ws(Twb)  
        
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
        
        h_fg = self.__properties.get_enthalpy_vaporization()                    # Latent Heat of Vaporization (Btu/lbm)        
        Cp_da = self.__properties.get_const_Cp_dry_air()                        # specific heat of dry air
        Cp_w = self.__properties.get_const_Cp_water()                           #  specific heat of water 
        Cp_v = self.__properties.get_const_Cp_vapor()                           # specific heat of vapor 
        Ws_wb = self.sat_humidity_ratio                                         # saturated humidity ratio function (pointer) Ws(Twb)  
        
        a = ( ( h_fg - ( Cp_w - Cp_v ) * T_wb ) * Ws_wb( T_wb ) - Cp_da * ( T_db - T_wb ) )
        b = ( h_fg + Cp_v * T_db - Cp_w * T_wb )
        
        return ( a / b )
    
    def grains_moisture( self ):
        pass
    
    def specific_humidity( self ):
        pass
    
    def absolute_humidity( self ):
        pass
    
    def degree_saturation( self, T_db, RH ):
        '''
        '''
        
        W = self.humidity_ratio( T_db, RH )
        Ws = self.sat_humidity_ratio( T_db )
        
        return ( W / Ws * 100)
    
    def specific_volume( self, T_db, RH ):
        '''
        '''
        T_db_R = self.convert.F_to_R( T_db )
        W = self.humidity_ratio(T_db, RH)
        P = self.__abs_pressure
        
        return 0.3704 * T_db_R * ( 1 + 1.6078 * W ) / P
    
    def density( self, T_db, RH ):
        
        return 1 / self.specific_volume(T_db, RH)
    
    def sensible_enthalpy( self ):
        pass
    
    def latent_enthalpy( self ):
        pass
    
    def total_enthalpy( self ):
        pass
    
    def get_pressure(self):
        
        print( self.__abs_pressure )


def main():
    pass

    x = PsychroCalcs( elevation = 5280 )
    
    print( "T_wb_reg: %.4f " %x.T_wb_regression( 50, 0.55 ) )
    print( "T_wb_iter: %.4f " %x.T_wb_iter( 50, 0.55 ) )
    print( "T_dp: %.4f " %x.T_dp( 50, 0.55 ) )
    print( "Hum Ratio: %.6f " %x.humidity_ratio( 50, 0.55 ) )
    print( "Hum Rat f(): %.6f " %x.W0( 50, 42.5958 ) )
    print( "Hum Rat f(): %.6f " %x.W1( 50, 42.5958 ) )
    print( "Density %.4f " %x.density( 50, 0.55 ) )
    
    x.get_pressure()


if __name__ == "__main__": 
    main()


