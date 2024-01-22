import sys

sys.path.append('../')

#from src.Pychrometrics import Pychrometrics
from PyChrometrics.src.Psychrometrics import Pychrometrics
from PyChrometrics.src.ThermodynamicProperties import MoistAirProperties

from PyChrometrics.src.PintSetup import Quantity

def main():

    qty = Quantity
    ma_props = MoistAirProperties()

    # elevation = qty(5280, "ft")
    # x = Pychrometrics( input_val = elevation, input_type = "elevation" )
    pressure = qty(12.5, "psi")
    x = Pychrometrics(input_val = pressure, input_type = "pressure")

    rh = qty(0.70, '%')
    Tdb = qty(68, 'degF')
    Twb = x.T_wb_iter(Tdb, rh)

    print(f'Twb: {Twb : 0.4f}')
    print(f'Tdp: {x.T_dp(Tdb, rh) : 0.4f}')
    print(f'W0: {x.W0(Tdb, Twb) : 0.4f}')

    
    # print( "T_wb_reg: %.4f " %x.T_wb_regression( Tdb, rh ) )
    # print( "T_wb_iter: %.4f " %x.T_wb_iter( Tdb, rh ) )
    # print( "T_dp: %.4f " %x.T_dp( Tdb, rh ) )
    # print( "Hum Ratio: %.6f " %x.W( Tdb, rh ) )
    # print( "Hum Rat f(): %.6f " %x.W0( Tdb, Twb ) )
    # print( "Hum Rat f(): %.6f " %x.W1( Tdb, Twb ) )
    # print( "Density %.4f " %x.density( Tdb, rh ) )
    # print( "Enthalpy %.4f "%x.enthalpy( Tdb, rh) )
    # print( "RH %.4f "%x.RH(Tdb, Tdp) )

    # y = x.W(65.3, 0.61) - x.W(59.2, 0.690)
    
    # print( "Delta W {:4f}".format(y) )
    
if __name__ == "__main__": 
    main()