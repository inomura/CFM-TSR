'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'C'

dHA = 17056.392158200473       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 67.00757392524397    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity
            
dH = 5347.447788538372      #Width of T in P      #Extremely sensitive

n = 614.9059778104743      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive
               
c = 6.224387566440748     #"Normal" sensitivity (linear to Mu)

Cd = 0.0039040554431937695

mDay = 0.164931474195987

Vmax = 3255.875133 #(um3)

O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
