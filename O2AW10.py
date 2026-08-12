'''
Created on Jun 8, 2018
This one has maintenance
Look at Kei210~111(123)-114(126) Kei211~10
@author: Keisuke
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *
from TemperatureViscosity import *   #in the source folder "Functions"
from DissolvedO2Saturation import *  #in the source folder "Functions"
from af003_energy_calculation import *
import random


#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Sub funcations
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

# Counting non-nan
def count_nonnan(a):
    return size(a) - count_nonzero(isnan(a))

# Getting average with nonnan values
def ave_nonnan(a):
    return nansum(a)/count_nonnan(a)

# giving perturbation 
def step1(A,i): #For smaller of the change
    if i > 0:
        while True:
            A0 = A + A/R*random.uniform(-1,1)
            if A0 > 0:
                break
    elif i == 0:
        A0 = A
    return A0

# giving perturbation 
def step2(A,i): #For bigger of the change
    if i > 0:
        while True:
            A0 = A + A/R2*random.uniform(-1,1)
            if A0 > 0:
                break
    elif i == 0:
        A0 = A
    return A0

# This is to obtain V based on best values for growth rate fit
def forVinO2(T,Kelvin,dHAbest,dCpBest,dHbest,nBest,cBest):    #Here we use O2 intead of CO2 for V computation (
                    #mainly for heterotrophic bacteria) modification largelly based on c106 03 01
    TK = T + Kelvin
    O2 = zeros(size(T))
    mu = zeros(size(T))
    i = 0
    for temp in T:  
        O2[i] = DissolvedO2Saturation(temp)/1000   #(mol O2 m-3) 
        mu[i] = TemperatureViscosity(temp + Kelvin)
        i = i + 1
    Do2_25=2.12e-9 #Diffusion coefficient of O2 in the water at 25C (m2/s) 
    mu25 = 890.3
    Tfactor = TK/(25+Kelvin)/mu*mu25  #Lerned from 633 00 07
    Do2 = Do2_25*Tfactor 
    pV = 1/(1+exp((-nBest)*(dHbest-18.1*TK+dCpBest*((TK-373.6)-TK*log(TK/385.2)))/(8.314*TK)))
    MuV = cBest*TK*exp(-dHAbest/(8.314*TK))*pV    #(d-1)
    MuV = MuV/86400     #(s-1)
    return O2,Do2,MuV

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Main Funcation
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
def O2AW(Alphabet,RT1,RT2,dHA,dCp,dH,n,c,Cd,mDay,Vmax):
    Vmax = Vmax * 1e-18 #Converting from (um3) to (m3)
    Vdata0 = genfromtxt('Weisse1998\\'+Alphabet+'-V.txt',delimiter=',')
    Mudata0 = genfromtxt('Weisse1998\\'+Alphabet+'-Mu.csv',delimiter=',')
    
    Kelvin = 273.15
    Vdata = Vdata0[:,1]*1e-18
    Mudata = Mudata0[:,1]
    
    TV = Vdata0[:,0]    #T for volume data
    TVK = TV + Kelvin
    TMu = Mudata0[:,0]   #T for Mu data
    TMuK = TMu + Kelvin
    
    MuAv=nansum(Mudata)/count_nonnan(Mudata)
    Vav = nansum(Vdata)/count_nonnan(Vdata)
    
    T = arange(0,35+1,0.1)
    TK = T + Kelvin  #Temperature in Kelvin
    p = 1/(1+exp((-n)*(dH-18.1*TK+dCp*((T-373.6)-TK*log(TK/385.2)))/(8.314*TK)))
    Mu = c*TK*exp(-dHA/(8.314*TK))*p
    Mu[Mu<0] = 0
    
    #AW method#===============Leraning from a801 20 00=========
    global R,R2
    global T0best,Abest,Cbest,Mumaxbest,tau2best,Bbest
    
    R = 300
    R2 = 10000
    
    #Rpeating part --------------------------------------------
    
    RpeatNum = arange(RT1)
    o = zeros(size(RpeatNum))

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # AW method for Mu
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    for i in RpeatNum:
        
        dHA0 = step1(dHA,i)
        dCp0 = step1(dCp,i)
        dH0 = step2(dH,i)
        n0 = step1(n,i)
        c0 = step1(c,i)
        
        p0 = 1/(1+exp((-n0)*(dH0-18.1*TMuK+dCp0*((TMuK-373.6)-TMuK*log(TMuK/385.2)))/(8.314*TMuK)))
        Mu0 = c*TMuK*exp(-dHA0/(8.314*TMuK))*p0
        
        X2 = (Mu0 - Mudata)**2/(2*MuAv**2)
    
        X2 = nansum(X2)/count_nonnan(X2)
        
        P_new=exp(-X2)
        
        
        if i == 0:
            dHA = dHA0
            dCp = dCp0
            dH = dH0
            n = n0
            c = c0
            
            dHAbest = dHA
            dCpBest = dCp
            dHbest = dH
            nBest = n
            cBest = c
            
            P = P_new
            Pbest = P_new
            
        
        else:
            Pratio = P_new/P
            r01 = random.uniform(0,1)
          #  print(P_new)
            
            if r01 < Pratio:
                P = P_new
                dHA = dHA0
                dCp = dCp0
                dH = dH0
                n = n0
                c = c0
                
                if P > Pbest:
                    Pbest = P
                    dHAbest = dHA
                    dCpBest = dCp
                    dHbest = dH
                    nBest = n
                    cBest = c
    
                    print('Pbest = ',Pbest)
    
    print('final Pbest = ',Pbest)    
    PbestMu = Pbest
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Mu for plotting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    p = 1/(1+exp((-nBest)*(dHbest-18.1*TK+dCpBest*((TK-373.6)-TK*log(TK/385.2)))/(8.314*TK)))
    Mu = cBest*TK*exp(-dHAbest/(8.314*TK))*p
    
    Mu[Mu<0] = 0
    #R2Mu computation##########################
    
    pR2 = 1/(1+exp((-nBest)*(dHbest-18.1*TMuK+dCpBest*((TMuK-373.6)-TMuK*log(TMuK/385.2)))/(8.314*TMuK)))
    MuR2 = cBest*TMuK*exp(-dHAbest/(8.314*TMuK))*pR2
    
    R2Mu = 1 - nansum((Mudata - MuR2)**2)/nansum((Mudata - MuAv)**2)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    O2,Do2,MuV = forVinO2(TV,Kelvin,dHAbest,dCpBest,dHbest,nBest,cBest)
    
    a0 = 0.216
    b = 0.939
    a = a0*10**(18*b)/(12*10**12)  #Unit converted for Qc (mol) - V (m3) conversion from Menden-Deuer 2000)
    
    RepeatNum2 = arange(RT2)
    Yo2C = McCarty().Yo2C
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # AW method for V
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    m = mDay/86400
    
    for i in RepeatNum2:
        
        Cd0 = step1(Cd,i)
        mDay0 = step1(mDay,i)
        
        m = mDay/86400
        
        V0 = ((48*pi**2)**(1/3)/a*Cd0*Do2*O2/((MuV + m)*Yo2C))**(1/(b-1/3))
        V0[V0>Vmax] = Vmax

        #Key part of AW method
        X2 = (V0 - Vdata)**2/(2*Vav**2)
        X2 = nansum(X2)/count_nonnan(X2)
        
        P_new=exp(-X2)
        ######################
        
        if i == 0:
            Cd = Cd0
            mDay = mDay0
            Cdbest = Cd
            mDaybest = mDay
            
            P = P_new
            Pbest = P_new
        
        else:
            Pratio = P_new/P
            r01 = random.uniform(0,1)
            
            if r01 < Pratio:
                Cd = Cd0
                mDay = mDay0
                
                P = P_new
                
                if P > Pbest:
                    Pbest = P
                    Cdbest = Cd
                    mDaybest = mDay
                    
                    print('Pbest2 = ',Pbest)
                    
    print('final Pbest2 = ', Pbest)
    PbestV = Pbest
    print('Cdbest = ',Cdbest)
    print('mDaybest = ',mDaybest)

    print('dHAbest = ',dHAbest)
    print('dCpBest = ',dCpBest)
    print('dHbest = ',dHbest)
    print('nBest = ',nBest)
    print('cBest = ',cBest)
    
    m = mDaybest/86400 #this is important!!!
  #  print('Vav = ',Vav*1e18)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # For V plotting
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    O2,Do2,MuV = forVinO2(T,Kelvin,dHAbest,dCpBest,dHbest,nBest,cBest)
    V = ((48*pi**2)**(1/3)/a*Cdbest*Do2*O2/((MuV + m)*Yo2C))**(1/(b-1/3))
    V[V==0]=nan #Here due to mu going to infinity (due to interpolation equation of ln), Do2 goes to zero thus V goes to zero. Here I call it nan so that it does not appear in the plot
    V[V>Vmax] = Vmax
    #R2V computation###########################
    O2,Do2,MuV = forVinO2(TV,Kelvin,dHAbest,dCpBest,dHbest,nBest,cBest)
    Vr2 = ((48*pi**2)**(1/3)/a*Cdbest*Do2*O2/(MuV + m)*Yo2C)**(1/(b-1/3))
    Vr2[Vr2>Vmax] = Vmax
    
    R2V = 1 - nansum((Vdata - Vr2)**2)/nansum((Vdata - Vav)**2)  #see https://bellcurve.jp/statistics/course/9706.html
    ###########################################
    
    ###################OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    # Plotting
    ###################OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Xlabel = '$\mathit{T}$ (C$^\circ$)'
    Xmin = 0
    Xmax = max(T)+1e-10
    
    Folder = 'T-V\\Weisse1998\\03MandVmax'
    DPI = 300
    
    figure(1)
    Ymax=max(Mudata*1.2)
    plot(T,Mu,color='green',linewidth=5)
    plot(TMu,Mudata,'o',color='green',mec='black',markersize=20)
    xlabel(Xlabel)
    ylabel('$\mathit{\mu}$ (d$^{-1}$)')
    xlim(Xmin,Xmax)
    ylim(ymax=Ymax)
  #  text(Xmax*0.05,Ymax*0.95,'$P$ = '+str(round(PbestMu,2)),verticalalignment='top')
    text(Xmax*(-0.25),Ymax*1.05,Alphabet,verticalalignment='top',fontsize = 38)
    Savefig3(Folder,Alphabet+'-Mu',DPI)
    
    figure(2)
    Ymax=max(Vdata*1e18)*1.3
    plot(T,V*1e18,color='#663300',linewidth=5)
    plot(TV,Vdata*1e18,'o',color='#663300',mec='black',markersize=20)
    xlim(Xmin,Xmax)
    ylim(ymin=0,ymax=Ymax)
    xlabel(Xlabel)
    ylabel('$\mathit{V}$ $\mu$m$^3$')
   # text(Xmax*0.05,Ymax*0.05,'$Cd$ = '+str(round(Cdbest,3)))
   # text(Xmax*0.95,Ymax*0.05,'$P$ = '+str(round(PbestV,2)),horizontalalignment='right',)
    text(Xmax*(-0.25),Ymax*1.05,Alphabet,verticalalignment='top',fontsize = 38)
    Savefig3(Folder,Alphabet+'-V',DPI)
    
    show()
