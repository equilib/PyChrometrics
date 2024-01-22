'''

MoistAirProperties - Imperial (English Units)

Library of typical functions that contain properties of moist air

** Dependencies **
1. Conversions.py - library of functions to convert between units

init development - 30/09/2017

@author: Gregory Evans
@version: 1.0.0

'''

from PyChrometrics.src.Conversions import UnitConversions
#from src.Conversions import UnitConversions

from PyChrometrics.src.PintSetup import Quantity
#from src.PintSetup import Quantity

qty = Quantity

class MoistAirProperties():
    
    #__Rda = 53.352                      # ft-lbf/(lb - R)
    __Rda = qty(53.352, 'ft * lbf / (lb * degR)')
    #__Rv = 85.778                       # ft-lbf/(lb - R)
    __Rv = qty(85.778, 'ft * lbf / (lb * degR)')
    #__cp_da = 0.240                     # Btu/(lbm - R) - specific heat capacity of dry air
    __cp_da = qty(0.240, 'Btu / (lb * degR)')
    #__cp_w = 0.999                      # Btu/(lbm - R) - specific heat capacity of water
    __cp_w = qty(0.999, 'Btu / (lb * degR)')
    #__cp_v = 0.451                      # Btu/(lbm - R) - specific heat capacity of vapor
    __cp_v = qty(0.451, 'Btu / (lb * degR)')
    #  __h_fg = 970.33                     # enthalpy of vaporization: Btu/lbm - Per Engineering Thermodynamics Steam Tables - Boles
    #__h_fg = 1060.9                     # enthalpy from Engineering Toolbox https://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html   *** WHERE DID THIS VALUE COME FROM? ***
    __h_fg = qty(1069.9, 'Btu / lb')
    #__molec_mass_ratio = 0.62198        # molecular mass ratio: Mw/Mda (ratio of molecular masses for vapor & dry air)
    __molec_mass_ratio = qty(0.62198)

    
    __GRAINS_PER_LBM = qty(7000,'grains')               # 7000 grains of moisture per lbm of air
    __FREEZING_POINT = qty(32.0, 'degF')                # deg F - freezing point of pure water
    #__STD_ATM__PRESSURE = qty(14.696, 'psi')

    def __init__(self):
        pass
    
    def R_da( self ) -> qty:
        return self.__Rda
    
    def R_v( self ) -> qty:
        return self.__Rv
    
    def Cp_da_const( self ) -> qty:
        return self.__cp_da
    
    def Cp_w_const( self ) -> qty:
        return self.__cp_w
    
    def Cp_v_const( self ) -> qty:
        return self.__cp_v
    
    def h_fg( self ) -> qty:
        return self.__h_fg
    
    def MMR( self ) -> qty:
        '''
        Return molecular mass ratio
        '''
        return self.__molec_mass_ratio
    
    def FREEZING_POINT( self ) -> qty:
        return self.__FREEZING_POINT

    def GRAINS_PER_LBM( self ) -> qty:
        return self.__GRAINS_PER_LBM
    

    def Cp_water(self, 
                 T_wb : qty ) -> qty:
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
        T_wb = T_wb.to('degC')
        
        # convert from J/kg-C to kJ/kg-C
        Cp_w = ( 0.0265 * T_wb.m ** 2 - 1.7688 * T_wb.m + 4205.6 ) * 1e-3
        Cp_w = qty(Cp_w, 'kJ / (kg * K)')
        
        return Cp_w.to('Btu / (lb * degF)')
    
    
    def Cp_vapor( self, 
                  T_db : qty ) -> qty:
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
        T_db = T_db.to('degC')

        # convert from J/kg-C to kJ/kg-C
        Cp_vapor = ( 0.0016 * T_db.m ** 2 + 0.1546 * T_db.m + 1858.7 ) * 1e-3
        Cp_vapor = qty(Cp_vapor, 'kJ / (kg * K)')
        
        return Cp_vapor.to('Btu / (lb * degF)')
    
    
    def Cp_dry_air( self, 
                    T_db : qty, 
                    T_wb : qty ) -> qty:
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
        T_wb = T_wb.to('degC')
        T_db = T_db.to('degC')
        
        # convert from J/kg-C to kJ/kg-C
        Cp_dry_air = ( 0.0667 * ( T_db.m + T_wb.m ) / 2 + 1005 ) * 1e-3
        Cp_dry_air = qty(Cp_dry_air, 'kJ / (kg * degC)')
        
        return Cp_dry_air.to('Btu / (lb * degF)')



