from thermop import thermoppt 
from thermop import adsTS_E 
from thermop import PES 
from pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from kin_temp0 import kineticparameter, dydt,solve_odes,solve_odes2,solve_odes3,solve_odes4
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
from thermop import specie_index
from math import log10 as log





# GOOD FOR OPTIMIZATION /// FOR ESTIMATING EFFECT OF  /// TIME, TEMP, CAT_MASS, EXIT_RATE, FEED_CONC /// ON THE KINETICS RATE FOR DIFF CASES

def time_effect(temperature_selected, moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  
    # for generic time range of 24 hrs
    title0=title_of_molecule  # name of the molecule studied
    TIME = 24 # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = temperature_selected # for each molecule study
    ppt = moleculer_property_function # ppt1, ppt2,
    
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    G_m, int_id, S_m, H_m = ppt(temp)
    tt0, y = solve_odes2(temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x)
    i=0
    tt2=len(tt0)
    tt=[]
    for i in range(tt2):
        tt.append(tt0[i]/3600) # sec to hr # that convert it from seconds into hours
    
    rr = y[:,6]
    hh = y[:,7]
    pp = y[:,8] 
    xxcat = y[:,0]/S_x 
    rrx = y[:,1]/S_x
    uux = y[:,2]/S_x
    hhx = y[:,3]/S_x   
    vvx = y[:,4]/S_x
    ppx = y[:,5]/S_x    
    print(rr,hh,pp,tt)
    max_p = max(pp)  # Find the maximum y value
    print('max_p = ', max_p)
    #max_T = tt[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    #print (max_T, round(max_p,5))
  
    plt.figure(2)       # <=========
    plt.title(title0)        # <=========    
    plt.plot(tt,rr, label='$R(T)$')
    plt.plot(tt,pp, label='$P(T)$')
    plt.plot(tt,hh, label='$H_2(T)$')
    plt.xlabel('Reaction Time in hours')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(3)       # <=========
    plt.title(title0)       # <=========
    plt.plot(tt,xxcat, label='$X(T)$')
    plt.plot(tt,rrx, label='$RX(T)$')
    plt.plot(tt,uux, label='$UX(T)$')
    plt.plot(tt,hhx, label='$HX(T)$')
    plt.plot(tt,vvx, label='$VX(T)$')
    plt.plot(tt,ppx, label='$PX(T)$')
    plt.xlabel('Reaction Time in hours')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return tt, rr, hh, pp, xxcat, rrx, uux, hhx, vvx, ppx    


def time_effect2(selected_time, temperature_selected, moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  
    #for specific time range 
    title0=title_of_molecule  # name of the molecule studied
    TIME = selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = temperature_selected # for each molecule study
    ppt = moleculer_property_function # ppt1, ppt2,
    
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    G_m, int_id, S_m, H_m = ppt(temp)
    tt0, y = solve_odes2(temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x)
    
    tt2=len(tt0)
    tt=[]
    i=0
    for i in range(tt2):
        tt.append(tt0[i]/1) # sec to hr # that convert it from seconds into hours
    
    rr = y[:,6]
    hh = y[:,7]
    pp = y[:,8] 
    xxcat = y[:,0]/S_x 
    rrx = y[:,1]/S_x
    uux = y[:,2]/S_x
    hhx = y[:,3]/S_x   
    vvx = y[:,4]/S_x
    ppx = y[:,5]/S_x    
    print(rr,hh,pp,tt)
    max_p = max(pp)  # Find the maximum y value
    print('max_p = ', max_p)
    #max_T = tt[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    #print (max_T, round(max_p,5))
  
    plt.figure(2)       # <=========
    plt.title(title0)        # <=========    
    plt.plot(tt,rr, label='$R(T)$')
    plt.plot(tt,pp, label='$P(T)$')
    plt.plot(tt,hh, label='$H_2(T)$')
    plt.xlabel('Reaction Time in sec')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(3)       # <=========
    plt.title(title0)       # <=========
    plt.plot(tt,xxcat, label='$X(T)$')
    plt.plot(tt,rrx, label='$RX(T)$')
    plt.plot(tt,uux, label='$UX(T)$')
    plt.plot(tt,hhx, label='$HX(T)$')
    plt.plot(tt,vvx, label='$VX(T)$')
    plt.plot(tt,ppx, label='$PX(T)$')
    plt.xlabel('Reaction Time in sec')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return tt, rr, hh, pp, xxcat, rrx, uux, hhx, vvx, ppx    

def time_effect3(selected_time, temperature_selected, moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  
    #for specific time range GOOD
    title0=title_of_molecule  # name of the molecule studied
    TIME = selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = temperature_selected # for each molecule study
    ppt = moleculer_property_function # ppt1, ppt2,
    
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    G_m, int_id, S_m, H_m = ppt(temp)
    tt0, y = solve_odes2(temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x)
    
    tt2=len(tt0)
    tt=[]
    xxcat=[]
    rrx=[]
    uux=[]
    hhx=[]
    vvx=[]
    ppx=[]
    rr=[]
    pp=[]
    hh=[]
    i=0
    
    for i in range(tt2):
        tt.append(tt0[i]/1) # sec to hr # that convert it from seconds into hours
    for i in range(tt2):
        xe = y[:,0]
        xxcat.append(float(xe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        rxe = y[:,1]
        rrx.append(float(rxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        uxe = y[:,2]
        uux.append(float(uxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        hxe = y[:,3]
        hhx.append(float(hxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        vxe = y[:,4]
        vvx.append(float(vxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        pxe = y[:,5]
        ppx.append(float(pxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        re = y[:,6]
        rr.append(float(re[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        he = y[:,7]
        hh.append(float(he[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        pe = y[:,8]
        pp.append(float(pe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
  
    print(rr,hh,pp,tt)
    max_p = max(pp)  # Find the maximum y value
    print('max_p = ', max_p)
    #max_T = tt[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    #print (max_T, round(max_p,5))
  
    plt.figure(2)       # <=========
    plt.title(title0)        # <=========    
    plt.plot(tt,rr, label='$R(T)$')
    plt.plot(tt,pp, label='$P(T)$')
    plt.plot(tt,hh, label='$H_2(T)$')
    plt.xlabel('Reaction Time in sec')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(3)       # <=========
    plt.title(title0)       # <=========
    plt.plot(tt,xxcat, label='$X(T)$')
    plt.plot(tt,rrx, label='$RX(T)$')
    plt.plot(tt,uux, label='$UX(T)$')
    plt.plot(tt,hhx, label='$HX(T)$')
    plt.plot(tt,vvx, label='$VX(T)$')
    plt.plot(tt,ppx, label='$PX(T)$')
    plt.xlabel('Reaction Time in sec')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return tt, rr, hh, pp, xxcat, rrx, uux, hhx, vvx, ppx    

#xxxxxxxxxxxxxxxx
title2='propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 

"""
# for generic time range of 24hr
tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect(525, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect(800, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect(625, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect(650, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect(700, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

# for specific time range peculiar to each molecules
tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect2(1.004,525, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect2(1.004,800, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect2(1.004,625, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect2(1.004,650, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect2(1.004,700, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

# for uniform time and change temp 400, 500 and 600K for each molecules
tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect2(1.004,600, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect2(1.004,600, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect2(1.004,600, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect2(1.004,600, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect2(1.004,600, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])
#OR 
tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect3(1.004,600, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect3(1.004,600, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect3(1.004,600, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect3(1.004,600, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect3(1.004,600, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

max_p1 = max(pp1)  # Find the maximum y value
max_p2 = max(pp2)  # Find the maximum y value
max_p3 = max(pp3)  # Find the maximum y value
max_p4 = max(pp4)  # Find the maximum y value
max_p5 = max(pp5)  # Find the maximum y value

print('max_p1=',max_p1)
print('max_p2=',max_p2)
print('max_p3=',max_p3)
print('max_p4=',max_p4)
print('max_p5=',max_p5)

max_T1 = np.interp(float(max_p1), pp1, tt1)
max_T2 = np.interp(float(max_p2), pp2, tt2)
max_T3 = np.interp(float(max_p3), pp3, tt3)
max_T4 = np.interp(float(max_p4), pp4, tt4)
max_T5 = np.interp(float(max_p5), pp5, tt5)

max_T1 = np.matrix(tt1)[np.matrix(pp1).index(max_p1)]  # Find the x value corresponding to the maximum y value
max_T2 = np.matrix(tt2)[np.matrix(pp2).index(max_p2)]  # Find the x value corresponding to the maximum y value
max_T3 = np.matrix(tt3)[np.matrix(pp3).index(max_p3)]  # Find the x value corresponding to the maximum y value
max_T4 = np.matrix(tt4)[np.matrix(pp4).index(max_p4)]  # Find the x value corresponding to the maximum y value
max_T5 = np.matrix(tt5)[np.matrix(pp5).index(max_p5)]  # Find the x value corresponding to the maximum y value

max_T1 = tt1[pp1.index(max_p1)]  # Find the x value corresponding to the maximum y value
max_T2 = tt2[pp2.index(max_p2)]  # Find the x value corresponding to the maximum y value
max_T3 = tt3[pp3.index(max_p3)]  # Find the x value corresponding to the maximum y value
max_T4 = tt4[pp4.index(max_p4)]  # Find the x value corresponding to the maximum y value
max_T5 = tt5[pp5.index(max_p5)]  # Find the x value corresponding to the maximum y value

print (title1,'time_max1=', max_T1, 'Pmax1=', round(max_p1,5))
print (title2,'time_max2=', max_T2, 'Pmax2=', round(max_p2,5))
print (title3,'time_max3=', max_T3, 'Pmax3=', round(max_p3,5))
print (title4,'time_max4=', max_T4, 'Pmax4=', round(max_p4,5))
print (title5,'time_max5=', max_T5, 'Pmax5=', round(max_p5,5))
"""
#xxxxxxxxxxxxxxxx



def time_effect4(mc0, C_R0, outflow_rate, selected_time, temperature_selected, moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  
    #for specific time range GOOD FOR OPTIMIZATION STUDY
    title0=title_of_molecule  # name of the molecule studied
    TIME = selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = temperature_selected # for each molecule study
    ppt = moleculer_property_function # ppt1, ppt2,
    
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    G_m, int_id, S_m, H_m = ppt(temp)
    tt0, y = solve_odes3(mc0,C_R0,outflow_rate,temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x)
    
    tt2=len(tt0)
    tt=[]
    xxcat=[]
    rrx=[]
    uux=[]
    hhx=[]
    vvx=[]
    ppx=[]
    rr=[]
    pp=[]
    hh=[]
    i=0
    
    for i in range(tt2):
        tt.append(tt0[i]/60) # sec to hr # that convert it from seconds into hours
    for i in range(tt2):
        xe = y[:,0]
        xxcat.append(float(xe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        rxe = y[:,1]
        rrx.append(float(rxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        uxe = y[:,2]
        uux.append(float(uxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        hxe = y[:,3]
        hhx.append(float(hxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        vxe = y[:,4]
        vvx.append(float(vxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        pxe = y[:,5]
        ppx.append(float(pxe[[i]])/S_x) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        re = y[:,6]
        rr.append(float(re[[i]])) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        he = y[:,7]
        hh.append(float(he[[i]])) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
    for i in range(tt2):
        pe = y[:,8]
        pp.append(float(pe[[i]])) # sec to hr # that convert it from seconds into hours    rr = y[:,6]
  
    plt.figure(2)       # <=========
    plt.title(title0)        # <=========    
    plt.plot(tt,rr, label='$R(T)$')
    plt.plot(tt,pp, label='$P(T)$')
    plt.plot(tt,hh, label='$H_2(T)$')
    plt.xlabel('Reaction Time in mins')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(3)       # <=========
    plt.title(title0)       # <=========
    plt.plot(tt,xxcat, label='$X(T)$')
    plt.plot(tt,rrx, label='$RX(T)$')
    plt.plot(tt,uux, label='$UX(T)$')
    plt.plot(tt,hhx, label='$HX(T)$')
    plt.plot(tt,vvx, label='$VX(T)$')
    plt.plot(tt,ppx, label='$PX(T)$')
    plt.xlabel('Reaction Time in mins')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return tt, rr, hh, pp, xxcat, rrx, uux, hhx, vvx, ppx    

title2='propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 

# the arguments are ====>   mc0,  C_R0,  outflow_rate,  Temp, time, dydt,  Gm, intID, specieID_list, S_x

"""
tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect4(25, 950, 5e-7, 1.004, 600, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect4(25, 950, 5e-7, 1.004, 600, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect4(25, 950, 5e-7, 1.004, 600, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect4(25, 950, 5e-7, 1.004, 600, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect4(25, 950, 5e-7, 1.004, 600, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])
print(rr1,hh1,pp1,tt1)
max_p1 = max(pp1)  # Find the maximum y value
print('max_p1 = ', max_p1)
"""

"""
#FOR OPTIMIZATION CONDITION

masss=30.4 # in mg
TIMMME=10 # in minutes
tempp=496.2 # in K
concc=1000  # in mol/m3
FLOWW=37.37  # in m3/sec

#XXXXXXXXXXXXXXXX 

timee=TIMMME*60 # in sec
rateflo=FLOWW*1e-7 # m3/sec

tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect4(masss, concc, rateflo, timee, tempp, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect4(masss, concc, rateflo, timee, tempp, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect4(masss, concc, rateflo, timee, tempp, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect4(masss, concc, rateflo, timee, tempp, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect4(masss, concc, rateflo, timee, tempp, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

#tt1,rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=time_effect4(30.4, 1000, 37.37e-7, 10*60, 496.2, ppt1, title1, ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2'])
#tt2,rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=time_effect4(masss, concc, rateflo, timee, tempp, ppt2, title2, ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2'])
#tt3,rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=time_effect4(masss, concc, rateflo, timee, tempp, ppt3, title3, ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2'])
#tt4,rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=time_effect4(masss, concc, rateflo, timee, tempp, ppt4, title4, ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2'])
#tt5,rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=time_effect4(masss, concc, rateflo, timee, tempp, ppt5, title5, ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2'])

ini_rr1=rr1[0]
ini_rr2=rr2[0]
ini_rr3=rr3[0]
ini_rr4=rr4[0]
ini_rr5=rr5[0]

max_p1 = max(pp1)  # Find the maximum y value
max_p2 = max(pp2)  # Find the maximum y value
max_p3 = max(pp3)  # Find the maximum y value
max_p4 = max(pp4)  # Find the maximum y value
max_p5 = max(pp5)  # Find the maximum y value

print('max_p1=',max_p1)
print('max_p2=',max_p2)
print('max_p3=',max_p3)
print('max_p4=',max_p4)
print('max_p5=',max_p5)

max_T1 = np.interp(float(max_p1), pp1, tt1)
max_T2 = np.interp(float(max_p2), pp2, tt2)
max_T3 = np.interp(float(max_p3), pp3, tt3)
max_T4 = np.interp(float(max_p4), pp4, tt4)
max_T5 = np.interp(float(max_p5), pp5, tt5)

print("------------------------------------------------------------------")
print (title1,'time_max1=', max_T1, 'Pmax1=', round(max_p1,5), 'conversion1=', round(max_p1/ini_rr1*100,2))
print (title2,'time_max2=', max_T2, 'Pmax2=', round(max_p2,5), 'conversion2=', round(max_p2/ini_rr2*100,2))
print (title3,'time_max3=', max_T3, 'Pmax3=', round(max_p3,5), 'conversion3=', round(max_p3/ini_rr3*100,2))
print (title4,'time_max4=', max_T4, 'Pmax4=', round(max_p4,5), 'conversion4=', round(max_p4/ini_rr4*100,2))
print (title5,'time_max5=', max_T5, 'Pmax5=', round(max_p5,5), 'conversion5=', round(max_p5/ini_rr5*100,2))
print("------------------------------------------------------------------")
print('masss=',masss,'TIMMME=',TIMMME, 'tempp=', tempp, 'concc=', concc, 'FLOWW=', FLOWW)

"""
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

