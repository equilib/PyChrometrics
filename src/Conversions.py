# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 07:27:46 2017

@author: gaevans
"""

class UnitConversions:    
    #--------------------------------------------------------------------------------------------------------------------------
    # TEMPERATURE CONVERSIONS
    # Rankine - returns the Thermodynamic temperature from deg F to deg R
    def F_to_R( self, T_farenheit ):
        return ( T_farenheit + 459.67 )
    
    # Kelvin - returns the Thermodynamic temperature from deg C to K
    def C_to_K( self, T_celsius ):
        return ( T_celsius + 273.15 )
    
    # F_to_C - returns the conversion of temperature from farenheit to celsius
    def F_to_C( self, T_farenheit ):
        return ( ( T_farenheit - 32.0 ) * 5.0 / 9.0 )
    
    # C_to_F - returns the conversion of temperature from celsius to farenheit
    def C_to_F( self, T_celsius ):
        return ( T_celsius * ( 9.0 / 5.0 ) + 32.0 )
    
    # K_to_C - returns the conversion of temperature from Kelvin to celsius
    def K_to_C( self, T_kelvin ):
        return ( T_kelvin - 273.15 )
    
    # R_to_F - returns the conversion of temperature from Rankine to farenheit
    def R_to_F( self, T_rankine ):
        return ( T_rankine - 459.67 )
    
    # Converts a given temperature to specified units
    
    def convert_temp( self, T, units_from = None, units_to = None ):
        
        if ( units_from != None and units_to != None ):
            if ( units_from.lower() == 'f' and units_to.lower() == 'c' ):
                return ( self.F_to_C( T ) )
            if ( units_from.lower() == 'c' and units_to.lower() == 'f' ):
                return ( self.C_to_F( T ) )
            if ( units_from.lower() == 'f' and units_to.lower() == 'k' ):
                return ( self.C_to_K( self.F_to_C( T ) ) )
            if ( units_from.lower() == 'c' and units_to.lower() == 'r' ):
                return ( self.F_to_R( self.C_to_F( T ) ) )
            if ( units_from.lower() == 'f' and units_to.lower() == 'r' ):
                return (self.F_to_R( T ) )
            if ( units_from.lower() == 'c' and units_to.lower() == 'k' ):
                return ( self.C_to_K( T ) )
        else: 
            print ( "Specify temperature units, none specified." )
    
        
    #--------------------------------------------------------------------------------------------------------------------------
    # PRESSURE CONVERSIONS
    
    # Pa_to_psia() - returns the conversion from Pascals to psi
    # @param P - pressure in pascals
    # @return pressure converted to psi
        
    def Pa_to_psia( self, P ):
        return ( P * 1.45038e-01 )
    
    def psia_to_Pa( self, P ):
        return ( P /1.45038e-01 )
    
    def kPa_to_psia( self, P ):
        return ( self.Pa_to_psia( P / 101.325 ) )
    
    def bar_to_psia( self, P ):
        return ( self, P * 1e3 )
    
    def mBar_to_psia( self, P ):
        return ( P * 0.145038 )
    
    
    #------------------------------------------------------------------------------------------------------------------------------
    # MISC Conversions
    def Btu_lbF_2_kJ_kgK( self, val ):
        return ( val * 4.1868 )           # kJ/kg-K
    
    def kJ_kgK_2_Btu_lbF( self, val ):
        return ( val * 0.238846 )         # Btu/lb-F
    
    def convert_RH( self, RH ):
        
        return( RH * 100 if RH < 1.0 and RH >= 0.0 else RH )
    