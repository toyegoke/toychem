from thermop import thermoppt 
from thermop import adsTS_E 
from thermop import PES 
from pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from kin_temp0 import kineticparameter, dydt,solve_odes
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
from thermop import specie_index

# GOOD FOR ESTIMATING EFFECT OF TEMP ON THE KINETICS RATE FOR DIFF CASES
"""
title2='propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 
TT = np.linspace(300, 1500, 2)  

def temp_effect1():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        G_m, int_id, S_m, H_m = ppt1(T)      # <============= 
        xx, xcat, rx, ux, hx, vx, px, r, h, p = solve_odes(T,0.25,dydt, G_m, int_id, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'],S_x)  # <===== 
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1])) 
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(1)       # <=========
    plt.title(title1)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(11)       # <=========
    plt.title(title1)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx    
"""


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


def temp_effect(temp, no_of_temp, moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  
    title0=title_of_molecule
    TT = np.linspace(350, temp, no_of_temp)    
    ppt = moleculer_property_function
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        G_m, int_id, S_m, H_m = ppt(T)      # <============= 
        xx, xcat, rx, ux, hx, vx, px, r, h, p = solve_odes(T,0.08333,dydt, G_m, int_id, Order_specieID_list,S_x)  # <===== 0.167 is 10min, 0.08333 is 5 min
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1])) 
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(2)       # <=========
    plt.title(title0)        # <=========    
    plt.plot(TT,rr, label='$R(T)$')
    plt.plot(TT,pp, label='$P(T)$')
    plt.plot(TT,hh, label='$H_2(T)$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(3)       # <=========
    plt.title(title0)       # <=========
    plt.plot(TT,xxcat, label='$X(T)$')
    plt.plot(TT,rrx, label='$RX(T)$')
    plt.plot(TT,uux, label='$UX(T)$')
    plt.plot(TT,hhx, label='$HX(T)$')
    plt.plot(TT,vvx, label='$VX(T)$')
    plt.plot(TT,ppx, label='$PX(T)$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return TT, rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx    

title2='propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 


TT1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=temp_effect(525, 10, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
TT2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=temp_effect(800, 10, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
TT3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=temp_effect(625, 10, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
TT4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=temp_effect(650, 10, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
TT5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=temp_effect(700, 10, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

TT1
TT2
TT3
TT4
TT5

max_p1 = max(pp1)  # Find the maximum y value
max_p2 = max(pp2)  # Find the maximum y value
max_p3 = max(pp3)  # Find the maximum y value
max_p4 = max(pp4)  # Find the maximum y value
max_p5 = max(pp5)  # Find the maximum y value


max_T1 = TT1[pp1.index(max_p1)]  # Find the x value corresponding to the maximum y value
max_T2 = TT2[pp2.index(max_p2)]  # Find the x value corresponding to the maximum y value
max_T3 = TT3[pp3.index(max_p3)]  # Find the x value corresponding to the maximum y value
max_T4 = TT4[pp4.index(max_p4)]  # Find the x value corresponding to the maximum y value
max_T5 = TT5[pp5.index(max_p5)]  # Find the x value corresponding to the maximum y value


print (title1,'Tmax1=', max_T1, 'Pmax1=', round(max_p1,5))
print (title2,'Tmax2=', max_T2, 'Pmax2=', round(max_p2,5))
print (title3,'Tmax3=', max_T3, 'Pmax3=', round(max_p3,5))
print (title4,'Tmax4=', max_T4, 'Pmax4=', round(max_p4,5))
print (title5,'Tmax5=', max_T5, 'Pmax5=', round(max_p5,5))





