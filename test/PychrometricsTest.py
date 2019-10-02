import sys

sys.path.append('../')

#from src.Pychrometrics import Pychrometrics
from PyChrometrics.src.Psychrometrics import Pychrometrics

def main():
    #pass

    x = Pychrometrics( elevation = 5600 )
    rh = 0.485
    Tdb = 68.0
    Twb = x.T_wb_iter(Tdb, rh)
    
    print( "T_wb_reg: %.4f " %x.T_wb_regression( Tdb, rh ) )
    print( "T_wb_iter: %.4f " %x.T_wb_iter( Tdb, rh ) )
    print( "T_dp: %.4f " %x.T_dp( Tdb, rh ) )
    print( "Hum Ratio: %.6f " %x.W( Tdb, rh ) )
    print( "Hum Rat f(): %.6f " %x.W0( Tdb, Twb ) )
    print( "Hum Rat f(): %.6f " %x.W1( Tdb, Twb ) )
    print( "Density %.4f " %x.density( Tdb, rh ) )

    # print( "Grains %.4f  " %x.grains_moisture(Tdb, rh) )

    # print( "Sens enthalpy %.4f  " %x.enthalpy_dry_air(Tdb, rh) )
    # print( "Late enthalpy %.4f  " %x.enthalpy_vapor(Tdb, rh) )
    
    # print( "Enthalpy: %.2f " %x.enthalpy(Tdb, rh) )

    # print( "Press atm: %.4f  " %x.P_atm_std() )

    # print( "Pw_sat: %.4f  " %x.P_w_sat(Tdb) )
    # print( "Pv_sat: %.4f  " %x.P_v_sat(Tdb) )

    e = x.P_v_partial(Tdb,rh)
    es = x.P_w_sat(Tdb)
    
    # print("RH: %.4f "%(e/es * 100) )

if __name__ == "__main__": 
    main()