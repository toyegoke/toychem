from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdatab
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5,ppt6
from BR_kin_temp import kineticparameter, dydtb,solve_odes46b
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import csv
from gn_thermop import specie_index
from math import log10 as log
from numpy import log as loge
from numpy import round
from numpy import polyfit

#=====================================================================================================
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS (BATCH)
#====================================================================================================

def RDS(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowrate):      # <=============  

    #for specific time range GOOD FOR RDS

    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K   
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime ## in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp ## in K # temperature_selected # for each molecule study
    mc0 = cmass  #=1e-3 ## in mg
    P_R0 = 1.01325e5*10 ## Pa 
    C_R0 = P_R0/R/temp ## mol/m3 
    outflow_rate = flowrate * (14.27) * 1e-7 ## in m3/sec
    v_R0 = (45*1e-6) ## in m3 (45 mL)
    n_R0 = v_R0*C_R0 ## in mol  

    ppt = moleculer_property_function ## ppt1, ppt2

    S_x_in = (18*1e-6) ## in mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) ## in umol/m^2 (equivalent of 0.44 ML from expt data)
     
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    # Co = float(p_o/R/298.15) # in mol/m^3
    
    #Catalyst properties
    dc = 29 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2    
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)     

    gat = [0,0,0,0,0,0,0,0,0,0,0];         bat = [1,1,1,1,1,1];         fat = bat;    rr0 = n_R0   

    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============    
        CCL = []        #RR = []
        RP = []
        RH = []
        T = temp
             
        xx1,y1= solve_odes46b(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydtb, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)

        xl1=len(xx1)
        x1=[]
        for i in range(xl1):
            x1.append(xx1[i]/3600) # sec to hr # that convert it from seconds into hours
        XX1=1-(y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5])
        rr1=y1[:,6]#*1e-6
        hh1=y1[:,7]#*1e-6
        pp1=y1[:,8]#*1e-6

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, y1[:,6]/y1[0,6], label='$n_R$')
        ax.semilogx(xx1, y1[:,7]/y1[0,6], label='$n_{H2,in}$')
        #ax.semilogx(xx1, (y1[:,8]-y1[:,7])/y1[0,6],  label='$n_{H2,out}$')
        ax.semilogx(xx1, y1[:,8]/y1[0,6], label='$n_P$')
        plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_1amount_plot_.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.semilogx(xx1, y1[:,0]/S_x, label='$\Theta_X$')
        ay.semilogx(xx1, y1[:,1]/S_x, label='$\Theta_{RX}$')
        ay.semilogx(xx1, y1[:,2]/S_x, label='$\Theta_{UX}$')
        ay.semilogx(xx1, y1[:,3]/S_x, label='$\Theta_{HX}$')
        ay.semilogx(xx1, y1[:,4]/S_x, label='$\Theta_{VX}$')
        ay.semilogx(xx1, y1[:,5]/S_x, label='$\Theta_{PX}$')
        plt.title(title+' Coverage Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
        ay.legend()  
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_1surface_plot_.png')

        # coverage and concentration definitions 1 (from back)
        n_R = y1[-1,6];  n_H2 = y1[-1,7]; n_P = y1[-1,8]
        Ttim = xx1[-1]; n_Ro = y1[0,6]

        # coverage and concentration definitions 2 (from back)
        n_Rb = y1[-2,6];  n_H2b = y1[-2,7]; n_Pb = y1[-2,8]
        Ttimb = xx1[-2]

        RRR=((n_R/n_Ro - n_Rb/n_Ro)/(Ttim - Ttimb)) 
        RPP=((n_P/n_Ro - n_Pb/n_Ro)/(Ttim - Ttimb)) 
        RHH=((n_H2/n_Ro - n_H2b/n_Ro)/(Ttim - Ttimb))

        glen=len(xx1)
        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_amount_data.cvs"
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerow(['Ri_mol/s'+str(T),'Pi_mol/s','Hi_mol/s','Time_s'])
            wr.writerow([ RRR*n_Ro, RPP*n_Ro, RHH*n_Ro ])
            wr.writerow(['Ri_1/s'+str(T),'Pi_1/s','Hi_1/s','Time_s'])
            wr.writerow([ RRR, RPP, RHH ])
            wr.writerow(['Ri/Ro'+str(T),'Hi/Ro ','Pi/Ro','Time'])
            i=0
            for i in range(glen):            
                wr.writerow([y1[i,6]/y1[0,6], y1[i,7]/y1[0,6], y1[i,8]/y1[0,6], xx1[i]])  
        return    
    return RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat)

#=====================================================================================================
#=====================================================================================================
# CALLING FOR SOLUTION TO BE RUN
#=====================================================================================================
#=====================================================================================================

title1= 'Ammonia Dehydrogenation over Ru'
engg1 = ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TS3','TSR','TSP','TSH2']
TIM=3.39
TIM2=TIM*1
pt = [ppt6]
titl = [title1]
specie_ID = [engg1]
#rtime = [10**(TIM)] # in sec 
rtime = [10**(5)] # in sec 
rtemp = [700+273.15]  # in Kelvin (K)
cmass = 50 # 1e-3 in mg
flowrate = 0 # batch will 0 while semibatch will be 1

ranP = len(titl)
i=0

#"""
for i in range(ranP):
        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowrate)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"""

i=1-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowrate)

"""

