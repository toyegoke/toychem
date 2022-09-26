from thermop import thermoppt 
from thermop import adsTS_E 
from thermop import PES, PESdata
from pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from kin_temp0 import kineticparameter, dydt,solve_odes,solve_odes2,solve_odes3,solve_odes4,solve_odes45
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import csv
from thermop import specie_index
from math import log10 as log
from numpy import log as loge
from numpy import round
from numpy import polyfit

#=====================================================================================================
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS 
#=====================================================================================================

def RDS(moleculer_property_function, title_of_molecule, Order_specieID_list):      # <=============  

    #for specific time range GOOD FOR RDS
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = (14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = 550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = 40  # in g
    C_R0 = 997.74 # mol/m3 
    outflow_rate = (14.27) *1e-7 # in m3/sec
    v_R0 = 45*1e-6 # in m3 (45 mL)
    n_R0 = v_R0*C_R0 # in mol  

    ppt = moleculer_property_function # ppt1, ppt2

    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    
    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    Co = float(p_o/R/298.15) # in mol/m^3
    
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2    
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)     

    gat = [0,0,0,0,0,0,0];         bat = [1,1,1,1,1,1];         fat = bat;    rr0 = n_R0   

    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============    
        CCL = []        #RR = []
        RP = []
        RH = []
        T = temp
             
        xx1,y1= solve_odes45(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)

        rxn_specie_id = Order_specieID_list     # list of selected specie IDs in a special order with string  elements
        x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
        vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
        tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13] 
        
        fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
        bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
        gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]
        
        energylist = G_m
        #Free Gibb (G) energy for the species
        Xg=energylist[specie_index(x,int_id)]
        Rg=energylist[specie_index(r,int_id)]
        RXg=energylist[specie_index(rx,int_id)]+gat1
        UXg=energylist[specie_index(ux,int_id)]+gat2
        HXg=energylist[specie_index(hx,int_id)]+gat3
        VXg=energylist[specie_index(vx,int_id)]+gat4
        PXg=energylist[specie_index(px,int_id)]+gat5
        H2g=energylist[specie_index(h2,int_id)]
        Pg=energylist[specie_index(p,int_id)]
        TS1g=energylist[specie_index(ts1,int_id)]+gat6
        TS2g=energylist[specie_index(ts2,int_id)]+gat7
        TSRg=energylist[specie_index(tsr,int_id)]
        TSPg=energylist[specie_index(tsp,int_id)]
        TSH2g=energylist[specie_index(tsh2,int_id)] 
        
        # calculate all reaction rate constants
        print('R+X<==>RX')
        #keqq1=kineticparameter.keq(Xg+Rg,RXg,T)*fat1/bat1
        kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T)*fat1 
        kb1=kineticparameter.ksu(TSRg,RXg,S_o,0,1,T)*bat1 
        #kb1=kf1/keqq1*Co
        
        print('RX+X<==>UX+HX')
        #keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
        kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)*fat2 
        kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)*bat2
        
        print('UX+X<==>VX+HX')
        #keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
        kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)*fat3 
        kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)*bat3
        
        print('VX<==>PX')
        #keqq3=kineticparameter.keq(PXg,VXg,T) 
        kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)*fat4 
        kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)*bat4  
        
        print('PX<==>P+X')
        #keqq4=kineticparameter.keq(Pg+Xg,PXg,T)*bat5/fat5 
        kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T)*bat5  
        kf4=kineticparameter.ksu(TSPg,PXg,S_o,0,1,T)*fat5  
        #kf4=kb4/keqq4*Co
        
        print('2HX<==>H2+2X')
        #keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)*bat6/fat6 
        kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,2,T)*bat6 
        kf5=kineticparameter.ksu(TSH2g,2*HXg,S_o,0,2,T)*fat6 
        #kf5=kb5/keqq5*Co  
        
        # coverage and concentration definitions
        X  = y1[-1,0]; RX = y1[-1,1]; UX = y1[-1,2] 
        VX = y1[-1,5]
        n_R = y1[-1,6];  n_H2 = y1[-1,7]; n_P = y1[-1,8]
	
	
        xl1=len(xx1)
        x1=[]
        for i in range(xl1):
            HX = y1[i,3]; PX = y1[i,4]
            x1.append(xx1[i]/3600) # sec to hr # that convert it from seconds into hours
            rP = SA*1e-6*kf4 * PX / n_R0  # in sec
            rH = SA*1e-6*kf5 * HX**2 /n_R0  # in sec
            #RR.append(rr)
            RP.append(rP)
            RH.append(rH)


        XX1=1-(y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5])
        rr1=y1[:,6]#*1e-6
        hh1=y1[:,7]#*1e-6
        pp1=y1[:,8]#*1e-6

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, RP, label='$n_R$')
        ax.semilogx(xx1, RH,marker='x', label='$n_{H2,out}$')
        plt.title(title+' desorpProduction Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_desorp_plot_.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, y1[:,6]/y1[0,6], label='$n_R$')
        ax.semilogx(xx1, y1[:,7]/y1[0,6],marker='x', label='$n_{H2,out}$')
        ax.semilogx(xx1, (y1[:,8]-y1[:,7])/y1[0,6], label='$n_{H2,in}$')
        ax.semilogx(xx1, y1[:,8]/y1[0,6],marker='x', label='$n_P$')
        plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_amount_plot_.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.semilogx(xx1, y1[:,0]/S_x,marker='x', label='$\Theta_X$')
        ay.semilogx(xx1, y1[:,1]/S_x, label='$\Theta_{RX}$')
        ay.semilogx(xx1, y1[:,2]/S_x,marker='x', label='$\Theta_{UX}$')
        ay.semilogx(xx1, y1[:,3]/S_x, label='$\Theta_{HX}$')
        ay.semilogx(xx1, y1[:,4]/S_x,marker='x', label='$\Theta_{VX}$')
        ay.semilogx(xx1, y1[:,5]/S_x, label='$\Theta_{PX}$')
        plt.title(title+' Coverage Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
        ay.legend()  
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_surface_plot_.png')

        print(title+' results')
        print('==================================')
        print('n(P)=',pp1[-1],'mol')
        print('n(H2)=',hh1[-1],'mol')
        print('n(R)=',rr1[-1],'mol')
        print('Y(P)=',pp1[-1]/rr0)
        print('Y(H2)=',hh1[-1]/rr0)
        print('X(R)=',(rr0-rr1[-1])/rr0)
        print('==================================')

        return        

    return RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat)

#=====================================================================================================
#=====================================================================================================
# CALLING FOR SOLUTION TO BE RUN
#=====================================================================================================
#=====================================================================================================

title1= '1-ol-oxacyclopentanol'
title2= 'Propanol'
title3= 'Cyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 

engg1 = ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2']
engg2 = ['X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2']
engg3 = ['X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2']
engg4 = ['X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2']
engg5 = ['X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2']

pt = [ppt1,ppt2,ppt3,ppt4,ppt5]
titl = [title1,title2,title3,title4,title5]
specie_ID = [engg1,engg2,engg3,engg4,engg5]
ranP = len(titl)
i=0

for i in range(ranP):

        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        RDS(ppt,title,engg)
        

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

