'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'E'

dHA = 45819.80968951132       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 64.39015385950839    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity
            
dH = 5340.510841530234     #Width of T in P      #Extremely sensitive

n = 1253.1638776820546      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive
           
c = 483885.75959961646     #"Normal" sensitivity (linear to Mu)

Cd = 0.0037340837679347883

mDay = 1.2782940517617212

Vmax = 2004.102433 #(um3)

O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
