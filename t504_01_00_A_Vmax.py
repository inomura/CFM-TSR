'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'A'

dHA = 13693.85578238447       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 65.5962719341015    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity

dH = 5339.656680288019      #Width of T in P      #Extremely sensitive

n = 1227.2628392947217      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive
         
c = 0.8981587432325384     #"Normal" sensitivity (linear to Mu)

Cd = 0.0017074984978005686

mDay = 0.11917418066041217

Vmax = 2112.5768


O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
