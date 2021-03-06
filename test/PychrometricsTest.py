import sys

sys.path.append('../')

#from src.Pychrometrics import Pychrometrics
from PyChrometrics.src.Psychrometrics import Pychrometrics

def main():

    x = Pychrometrics( elevation = 5600 )
    rh = 0.653
    Tdb = 64.3
    Twb = x.T_wb_iter(Tdb, rh)
    
    print( "T_wb_reg: %.4f " %x.T_wb_regression( Tdb, rh ) )
    print( "T_wb_iter: %.4f " %x.T_wb_iter( Tdb, rh ) )
    print( "T_dp: %.4f " %x.T_dp( Tdb, rh ) )
    print( "Hum Ratio: %.6f " %x.W( Tdb, rh ) )
    print( "Hum Rat f(): %.6f " %x.W0( Tdb, Twb ) )
    print( "Hum Rat f(): %.6f " %x.W1( Tdb, Twb ) )
    print( "Density %.4f " %x.density( Tdb, rh ) )

    y = x.W(65.3, 0.61) - x.W(59.2, 0.690)
    
    print( "Delta W {:4f}".format(y) )

if __name__ == "__main__": 
    main()