'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'F'

dHA = 22898.933016991577       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 65.38659810606394    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity
            
dH = 5341.061402754082     #Width of T in P      #Extremely sensitive

n = 1878.269315261679      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive
               
c = 42.14757924331265    #"No8rmal" sensitivity (linear to Mu)

Cd = 0.0025720464652130184

mDay = 0.4903354506541906

Vmax = 2147.592867 #(um3)

O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
