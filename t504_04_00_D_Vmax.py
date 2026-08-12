'''
Created on Jun 8, 2018
from t200 01 01
@author: Keisuke
'''


from O2AW10 import *

RT1=1#repeat times
RT2=1#repeat times

Alphabet = 'D'

dHA = 8302.227797538879       #Overall value and shape. Mu become high with low value. More flat with low value
                    #Shape is not so sensitive. However, the value is very sensitive -> should be compensated by "c"
                    
dCp = 65.31271219226387    #change the location of T range in P (high -> high T)
            #"Normal" sensitivity
            
dH = 5340.848400543525      #Width of T in P      #Extremely sensitive

n = 2017.271911243549      #Shape of the P (High -> more flat at max, low, value become lower but more wide range)
               #Somewhat sensitive

          
c = 0.08917155696255091    #"Normal" sensitivity (linear to Mu)

Cd = 0.001685486087145346

mDay = 7.153185103119407e-05

Vmax = 2267.78215 #(um3)

O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax)
