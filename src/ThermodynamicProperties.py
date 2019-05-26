'''

MoistAirProperties - Imperial (English Units)

Library of typical functions that contain properties of moist air

** Dependencies **
1. Conversions.py - library of functions to convert between units

init development - 30/09/2017

@author: Gregory Evans
@version: 1.0.0

'''

from Conversions import UnitConversions


class MoistAirProperties():
    
    convert = UnitConversions()
    
    __Rda = 53.352                      # ft-lbf/(lb - R)
    __Rv = 85.778                       # ft-lbf/(lb - R)
    __cp_da = 0.240                     # Btu/(lbm - R) - specific heat capacity of dry air
    __cp_w = 0.999                      # Btu/(lbm - R) - specific heat capacity of water
    __cp_v = 0.451                      # Btu/(lbm - R) - specific heat capacity of vapor
    #  __h_fg = 970.33                     # enthalpy of vaporization: Btu/lbm - Per Engineering Thermodynamics Steam Tables - Boles
    __h_fg = 1061                     # enthalpy from Engineering Toolbox https://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html   *** WHERE DID THIS VALUE COME FROM? ***
    __molec_mass_ratio = 0.62198        # molecular mass ratio: Mw/Mda (ratio of molecular masses for vapor & dry air)
    
    __GRAINS_PER_LBM = 7000             # 7000 grains of moisture per lbm of air
    __FREEZING_POINT = 32.0             # deg F - freezing point of pure water
    
    def get_R_dry_air( self ):
        return self.__Rda
    
    def get_R_vapor( self ):
        return self.__Rv
    
    def get_const_Cp_dry_air( self ):
        return self.__cp_da
    
    def get_const_Cp_water( self ):
        return self.__cp_w
    
    def get_const_Cp_vapor( self ):
        return self.__cp_v
    
    def get_enthalpy_vaporization( self ):
        return self.__h_fg
    
    def get_molec_mass_ratio( self ):
        return self.__molec_mass_ratio

    def h_fg( self ):
        return self.__h_fg
    
    def FREEZING_POINT( self ):
        return self.__FREEZING_POINT

    def GRAINS_PER_LBM( self ):
        return self.__GRAINS_PER_LBM
    

    def Cp_water( self, T_wb ):
        '''
        Cp_water - specific heat capacity of water
        
        Adapted from research listed in:
        A. M. Al-Ismaili, Modelling of a humidification dehumidification greenhouse in Oman.
        University, 2009, pp. 190 - 191.

        K. Raznjevic, Handbook of thermodynamic tables and charts
        Washington: Hemisphere Publishing Corp., 1975 pp. 67-68.
        
        @param - T_wb: wet bulb temperature (deg F)
        @return specific heat capacity of water (Btu/lbm-F)
        
        '''
        
        T_wb = self.convert.F_to_C( T_wb )
        
        # convert from J/kg-C to kJ/kg-C
        Cp_w = ( 0.0265 * T_wb ** 2 - 1.7688 * T_wb + 4205.6 ) * 1e-3
        
        return ( self.convert.kJ_kgK_2_Btu_lbF( Cp_w ) )
    
    
    def Cp_vapor( self, T_db ):
        '''
        Cp_vapor - specific heat capacity of water vapor
        
        Adapted from research listed in:
        A. M. Al-Ismaili, Modelling of a humidification dehumidification greenhouse in Oman.
        University, 2009, pp. 190 - 191.

        K. Raznjevic, Handbook of thermodynamic tables and charts
        Washington: Hemisphere Publishing Corp., 1975 pp. 67-68.
        
        @param - T_wb: wet bulb temperature (deg F)
        @return specific heat capacity of water vapor (Btu/lbm-F)
        
        '''
        
        T_db = self.convert.F_to_C( T_db )
        # convert from J/kg-C to kJ/kg-C
        Cp_vapor = ( 0.0016 * T_db ** 2 + 0.1546 * T_db + 1858.7 ) * 1e-3
        
        return ( self.convert.kJ_kgK_2_Btu_lbF( Cp_vapor ) )
    
    
    def Cp_dry_air( self, T_db, T_wb ):
        '''
        Cp_dry_air - specific heat capacity of dry air
        
        Adapted from research listed in:
        A. M. Al-Ismaili, Modelling of a humidification dehumidification greenhouse in Oman.
        University, 2009, pp. 190 - 191.

        K. Raznjevic, Handbook of thermodynamic tables and charts
        Washington: Hemisphere Publishing Corp., 1975 pp. 67-68.
        
        @param - T_db: dry bulb temperature (deg F)
        @param - T_wb: wet bulb temperature (deg F)
        @return specific heat capacity of water vapor (Btu/lbm-F)
        
        '''
        
        T_wb = self.convert.F_to_C( T_wb )
        T_db = self.convert.F_to_C( T_db )
        
        # convert from J/kg-C to kJ/kg-C
        Cp_dry_air = ( 0.0667 * ( T_db + T_wb ) / 2 + 1005 ) * 1e-3
        
        return ( self.convert.kJ_kgK_2_Btu_lbF( Cp_dry_air ) )
    
    
def main():
    
    props = MoistAirProperties()
    
    print( props.Cp_dry_air(50, 0.55) )
    
if __name__ == "__main__":
    main()


