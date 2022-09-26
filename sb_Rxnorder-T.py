from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdata
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from sb_kin_temp0 import kineticparameter, dydt,solve_odes,solve_odes2,solve_odes3,solve_odes4,solve_odes40
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

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# GOOD FOR RATE-DETERMING-STEP CALCULATION /// ON THE KINETICS RATE FOR DIFF CASES //  REACTION ORDER (N_ORDER-TEMP)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def EaApp(moleculer_property_function, title_of_molecule, Order_specieID_list, C_R0, reactiontime, reactiontemp, catmass, flowout):      # <=============  

    #for specific time range GOOD FOR REACTION ORDER (N_ORDER)
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = reactiontime #(14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = reactiontemp #550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = catmass #40  # in g
    #C_R0 = 997.74 # mol/m3 
    outflow_rate = flowout * (14.27) * 1e-7 # in m3/sec
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

    def EA_effect(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============  

        tt=[]; xxcat=[]; rrx=[]; uux=[]; hhx=[]; vvx=[]; ppx=[]
        rr1=[]; pp1=[]; hh1=[]; rr2=[]; pp2=[]; hh2=[]
        
        CCL = []        #RR = []
        RP = []
        RH = []
        
        CC = [C_R0-100,C_R0-80,C_R0-60,C_R0-40,C_R0-20,C_R0,C_R0+20,C_R0+40,C_R0+60,C_R0+80,C_R0+100]
        CC1 = [CC[0]/Co,CC[1]/Co,CC[2]/Co,CC[3]/Co,CC[4]/Co,CC[5]/Co,CC[6]/Co,CC[7]/Co,CC[8]/Co,CC[9]/Co,CC[10]/Co]        

        CCo=[]
        cclen=len(CC1)
        i=0

        for i in range(cclen):
            CCo.append(float(loge(CC1[i]))) # the nat log of the change in k i.e. ln(k_dff)

        tlens = len(CC)
        i = 0
        
        for i in range(tlens):
            
            T=temp
            C_R = CC[i]
            
            filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-Order-T/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
            filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-Order-T/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)
        
            G_m, int_id, S_m, H_m, fig, fig2 = ppt(T, filename_link_fomat, filename_link_fomat2)
        
            fig.savefig('/home/toyegoke/ENS_kinetics/Final-Order-T/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
            fig2.savefig('/home/toyegoke/ENS_kinetics/Final-Order-T/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)          
            
            tt0, y= solve_odes40(SA,n_R0,v_R0,C_R,outflow_rate,T,TIME,dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)
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
            X  = y[-1,0]; RX = y[-1,1]; UX = y[-1,2] 
            HX = y[-1,3]; PX = y[-1,4]; VX = y[-1,5]
            n_R = y[-1,6];  n_H2 = y[-1,7]; n_P = y[-1,8]
        
            rP = SA*1e-6*kf4 * PX / n_R0  # in sec
            rH = SA*1e-6*kf5 * HX**2 /n_R0  # in sec

            #RR.append(rr)

            RP.append(rP)
            RH.append(rH)
            CCL.append(C_R)
            
        return CCL, RH, RP, loge(RH), loge(RP), CCo       
                     #float(hh[-1]/rr[0]), float(pp[-1]/rr[0])    

    def gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list):
       
       title=title_of_molecule
       gat=[0,0,0,0,0,0,0]
       bat=[1,1,1,1,1,1]
       fat=bat
       ppt = moleculer_property_function # ppt1, ppt2,
       
       conc, RH, RP, lnRH, lnRP, conc0  = EA_effect(ppt, title0, Order_specieID_list, fat, bat, gat)
       slope1, intercept1 = np.polyfit(conc, lnRH, 1)
       slope2, intercept2 = np.polyfit(conc, lnRP, 1)
       slope11, intercept11 = np.polyfit(conc, RH, 1)
       slope22, intercept22 = np.polyfit(conc, RP, 1)       
       AppEa1 = C_R0*slope1      #C_R0 they are ORDERS not ACTIVATION ENERGIES
       AppEa2 = C_R0*slope2
       AppEa11 = C_R0*slope11
       AppEa22 = C_R0*slope22

       templen=len(conc) 
       i=0 
       lnRHfit=[]
       lnRPfit=[]
       RHfit=[]
       RPfit=[]
       
       for i in range(templen):
           aaa1 = slope1*float(conc[i]) + intercept1
           aaa2 = slope2*float(conc[i]) + intercept2
           aaa11 = slope11*float(conc[i]) + intercept11
           aaa22 = slope22*float(conc[i]) + intercept22
           lnRHfit.append(aaa1)
           lnRPfit.append(aaa2)
           RHfit.append(aaa11)
           RPfit.append(aaa22)

       return  AppEa1, AppEa2

    return gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list)


title1= '1-ol-oxacyclopentanol'
title2='propanol'
title3= 'cyclopentanol'
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

rtime = [3854.608589,666.2822896,3214.544089,2375.114529,1377.532526] # for Tm+25K +5%
#rtime = [30181.24358,36190.7871,5843905.391,18596.93624,2088453.228] # for Tm+25K +50%
rtemp = [502,460,418,502,539] # Tmin+25K
#rtime = [59.18063593,43.72650814,55.70480542,32.30799204,555.6449496] # for 550K +5%
#rtime = [463.3791336,798.9492025,901.7644281,303.3386229,702521.0872] # for 550K +50%
#rtemp = [550,550,550,550,550] # 550K

cmass = 1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

CONC = [197.74, 397.74, 597.74, 797.74, 997.74, 1097.74, 1197.74, 1297.74]
ranC = len(CONC)
ranP = len(titl)
i = 0

for i in range(ranP):

    ppt=pt[i]
    title=titl[i]
    engg=specie_ID[i]           #RDS=stepRDS[i]
    j=0
    Eaa1=[]
    Eaa2=[]

    for j in range(ranC):

        CR=CONC[j]
        AppEa1, AppEa2 = EaApp(ppt,title,engg,CR,rtime[i],rtemp[i],cmass,flowout)
        Eaa1.append(AppEa1)  
        Eaa2.append(AppEa2)          

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(CONC, Eaa1, label=title) #marker='o',
    plt.xlabel('$Conc / (mol/m^3)$')
    plt.ylabel('$Order / n_R$')
    plt.title(title+' Dehydrogenation')
    ax.legend()
    fig.savefig('/home/toyegoke/ENS_kinetics/Final-Order-T/'+title+'_plot_for_Order(T,H2).png')

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(CONC, Eaa2, label=title) #marker='o',
    plt.xlabel('$Conc / (mol/m^3)$')
    plt.ylabel('$Order / n_R$')
    plt.title(title+' Dehydrogenation')
    ax.legend()
    fig.savefig('/home/toyegoke/ENS_kinetics/Final-Order-T/'+title+'_plot_for_Order(T,P).png')
    
print('CALCULATION COMPLETED')
