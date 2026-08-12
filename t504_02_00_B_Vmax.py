'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'B'

dHA = 23112.305156157225       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 64.76545580487192    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity
            
dH = 5342.201489524019      #Width of T in P      #Extremely sensitive

n = 1505.3366110022473      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive
               
c = 31.990472118394575    #"Normal" sensitivity (linear to Mu)

Cd = 0.0016861919717323307

mDay = 0.514527506928996

Vmax = 1233.4052 #(um3)

O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
