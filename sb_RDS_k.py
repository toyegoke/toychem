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

# GOOD FOR RATE-DETERMING-STEP CALCULATION /// ON THE KINETICS RATE FOR DIFF CASES // APPARENT ACT_ENERGY (App.Ea)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def RDS(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout):      # <=============  

    #for specific time range GOOD FOR RDS (for RDS when rate is measured in terms of Rate_x)
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime # (14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp # 550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = cmass # 40  # in g
    C_R0 = 997.74 # mol/m3 
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
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-RDS-k/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-RDS-k/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)          

    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============    
        CCL = []        #RR = []
        RP = []; RPP = []
        RH = []; RHH = []; Rrr = []
        T = temp
            
        tt0, y= solve_odes46(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)
        """
        #========================================
        xx1=tt0  # time in sec
        
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, y[:,6]/y[0,6], label='$n_R$')
        ax.semilogx(xx1, y[:,7]/y[0,6], label='$n_{H2,in}$')
        ax.semilogx(xx1, (y[:,8]-y[:,7])/y[0,6],  label='$n_{H2,out}$')
        ax.semilogx(xx1, y[:,8]/y[0,6], label='$n_P$')
        plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_'+str(int(T))+'K_fat'+str(int(fat[0]))+'_1amount_plot_.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.semilogx(xx1, y[:,0]/S_x, label='$\Theta_X$')
        ay.semilogx(xx1, y[:,1]/S_x, label='$\Theta_{RX}$')
        ay.semilogx(xx1, y[:,2]/S_x, label='$\Theta_{UX}$')
        ay.semilogx(xx1, y[:,3]/S_x, label='$\Theta_{HX}$')
        ay.semilogx(xx1, y[:,4]/S_x, label='$\Theta_{VX}$')
        ay.semilogx(xx1, y[:,5]/S_x, label='$\Theta_{PX}$')
        plt.title(title+' Coverage Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
        ay.legend()  
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_'+str(int(T))+'K_fat'+str(int(fat[0]))+'_1surface_plot_.png')
        #========================================
        """
        rxn_specie_id = Order_specieID_list     # list of selected specie IDs in a special order with string  elements
        x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
        vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
        tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13] 

        fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
        bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
        gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]; gat8=gat[7]; gat9=gat[8]; gat0=gat[9]
    
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
        TSRg=energylist[specie_index(tsr,int_id)]+gat8
        TSPg=energylist[specie_index(tsp,int_id)]+gat9
        TSH2g=energylist[specie_index(tsh2,int_id)]+gat0         

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
        
        # coverage and concentration definitions 1 (from back)
        X  = y[-1,0]; RX = y[-1,1]; UX = y[-1,2] 
        HX = y[-1,3]; VX = y[-1,4]; PX = y[-1,5]
        n_R = y[-1,6];  n_H2 = y[-1,7]; n_P = y[-1,8]
        Ttim = tt0[-1]; n_Ro = y[0,6]

        # coverage and concentration definitions 2 (from back)
        Xb  = y[-2,0]; RXb = y[-2,1]; UXb = y[-2,2] 
        HXb = y[-2,3]; VXb = y[-2,4]; PXb = y[-2,5]
        n_Rb = y[-2,6];  n_H2b = y[-2,7]; n_Pb = y[-2,8]
        Ttimb = tt0[-2]
    
        rP = SA*1e-6*kf4 * PX / n_R0  # in sec
        rH = SA*1e-6*kf5 * HX**2 /n_R0  # in sec

        RP.append(rP)
        RH.append(rH)

        fatK = [kf1,kf2,kf3,kf6,kf4,kf5]
        print(fatK)    

        Rrr.append((n_R/n_Ro - n_Rb/n_Ro)/(Ttim - Ttimb)) 
        RPP.append((n_P/n_Ro - n_Pb/n_Ro)/(Ttim - Ttimb)) 
        RHH.append((n_H2/n_Ro - n_H2b/n_Ro)/(Ttim - Ttimb))

        #return RH[0], RP[0], loge(RH[0]), loge(RP[0]), fatK    
        return RHH[0], RPP[0], loge(RHH[0]), loge(RPP[0]), fatK     

    def gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list):
    
        # title0=titlemolecule
        dff = [0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9]

        gat = [0,0,0,0,0,0,0, 0,0,0]           #bat = [1,1,1,1,1,1]         #fat = bat 
        T = temp

        ppt = moleculer_property_function # ppt1, ppt2,

        rangedff=len(dff)
        
        rrr1=[]; ppp1=[]; hhh1=[]; rrr2=[]; ppp2=[]; hhh2=[] 
        rrr3=[]; ppp3=[]; hhh3=[]; rrr4=[]; ppp4=[]; hhh4=[] 
        rrr5=[]; ppp5=[]; hhh5=[]; rrr6=[]; ppp6=[]; hhh6=[]
        rrr7=[]; ppp7=[]; hhh7=[]
        i=0    
        lnG1=[]; lnG2=[]; lnG3=[]; lnG4=[]; lnG5=[]; lnG6=[]; lnG7=[]

        for i in range(rangedff):
        
            bat1=[dff[i],1,1,1,1,1] 
            bat2=[1,dff[i],1,1,1,1] 
            bat3=[1,1,dff[i],1,1,1] 
            bat4=[1,1,1,dff[i],1,1] 
            bat5=[1,1,1,1,dff[i],1] 
            bat6=[1,1,1,1,1,dff[i]] 
            
            fat1=bat1; fat2=bat2; fat3=bat3; fat4=bat4; fat5=bat5; fat6=bat6 
      
            # the arguments are ====>   mc0,  C_R0,  outflow_rate,  Temp, time, dydt,  Gm, intID, specieID_list, S_x, fat, bat, gat
            RH1, RP1, lnRH1, lnRP1, fatK1 = RDS_effect2(ppt, title0, Order_specieID_list, fat1, bat1, gat)
            RH2, RP2, lnRH2, lnRP2, fatK2 = RDS_effect2(ppt, title0, Order_specieID_list, fat2, bat2, gat)
            RH3, RP3, lnRH3, lnRP3, fatK3 = RDS_effect2(ppt, title0, Order_specieID_list, fat3, bat3, gat)
            RH4, RP4, lnRH4, lnRP4, fatK4 = RDS_effect2(ppt, title0, Order_specieID_list, fat4, bat4, gat)
            RH5, RP5, lnRH5, lnRP5, fatK5 = RDS_effect2(ppt, title0, Order_specieID_list, fat5, bat5, gat)
            RH6, RP6, lnRH6, lnRP6, fatK6 = RDS_effect2(ppt, title0, Order_specieID_list, fat6, bat6, gat)

            hhh1.append(lnRH1); ppp1.append(lnRP1) 
            hhh2.append(lnRH2); ppp2.append(lnRP2) 
            hhh3.append(lnRH3); ppp3.append(lnRP3) 
            hhh4.append(lnRH4); ppp4.append(lnRP4) 
            hhh5.append(lnRH5); ppp5.append(lnRP5) 
            hhh6.append(lnRH6); ppp6.append(lnRP6)
            
            lnG1.append(float(loge(fatK1[0]))) # the nat log of the change in G i.e. ln(G)
            lnG2.append(float(loge(fatK2[1]))) # the nat log of the change in G i.e. ln(G)
            lnG3.append(float(loge(fatK3[2]))) # the nat log of the change in G i.e. ln(G)
            lnG4.append(float(loge(fatK4[3]))) # the nat log of the change in G i.e. ln(G)
            lnG5.append(float(loge(fatK5[4]))) # the nat log of the change in G i.e. ln(G)
            lnG6.append(float(loge(fatK6[5]))) # the nat log of the change in G i.e. ln(G)

        lnG=[lnG1,lnG2,lnG3,lnG4,lnG5,lnG6]; gaps = dff 
        print(lnG)

        return hhh1,hhh2,hhh3,hhh4,hhh5,hhh6, ppp1,ppp2,ppp3,ppp4,ppp5,ppp6, lnG,gaps

    return gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list)
    
# solve for a and b // intercept & slope
    
def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))
    
    return a, b


title1= '1-ol-oxacyclopentanol'
title2= 'propanol'
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

"""
#rtime = [3854.608589,666.2822896,3214.544089,2375.114529,1377.532526] # for Tm+25K +5%
#rtime = [30181.24358,36190.7871,5843905.391,18596.93624,2088453.228] # for Tm+25K +50%
#rtemp = [502,460,418,502,539] # Tmin+25K
#rtime = [59.18063593,43.72650814,55.70480542,32.30799204,555.6449496] # for 550K +5%
#rtime = [463.3791336,798.9492025,901.7644281,303.3386229,702521.0872] # for 550K +50%
#rtemp = [550,550,550,550,550] # 550K
"""

cmass = 1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

#rtime=[3679.861095,658.6856963,2571.360722,2285.111635,1373.71884]# for Tm+25K +5%
#rtime=[13655.66496,3801.916238,17285.273,8486.623985,90089.66559]# for Tm+25K +20%
#rtemp = [502,460,418,502,539] # Tmin+25K

#rtime=[57.33126605,42.11744463,53.87014991,30.5556788,531.9821618]# for 550K +5%
#rtime=[214.6153297,247.5025416,294.6326122,127.2422,31021.55806]# for 550K +20%
rtemp = [550,550,550,550,550] # 550K
rtime = [10**(5),10**(5),10**(5),10**(5),10**(5)] #for 5%/20%

ranP = 1 # len(titl)
i=0

for i in range(ranP):

        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        hhh1,hhh2,hhh3,hhh4,hhh5,hhh6,ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,dff,gaps = RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout)
        print(ppp1, dff)
        
        slope1, intercept1 = np.polyfit(dff[0], ppp1, 1)
        slope2, intercept2 = np.polyfit(dff[1], ppp2, 1)
        slope3, intercept3 = np.polyfit(dff[2], ppp3, 1)
        slope4, intercept4 = np.polyfit(dff[3], ppp4, 1)
        slope5, intercept5 = np.polyfit(dff[4], ppp5, 1)
        slope6, intercept6 = np.polyfit(dff[5], ppp6, 1)

        yfit1 = [intercept1 + slope1 * float(xi) for xi in dff[0]]
        yfit2 = [intercept2 + slope2 * float(xi) for xi in dff[1]]
        yfit3 = [intercept3 + slope3 * float(xi) for xi in dff[2]]
        yfit4 = [intercept4 + slope4 * float(xi) for xi in dff[3]]
        yfit5 = [intercept5 + slope5 * float(xi) for xi in dff[4]]
        yfit6 = [intercept6 + slope6 * float(xi) for xi in dff[5]]

        slope11, intercept11 = np.polyfit(dff[0], hhh1, 1)
        slope22, intercept22 = np.polyfit(dff[1], hhh2, 1)
        slope33, intercept33 = np.polyfit(dff[2], hhh3, 1)
        slope44, intercept44 = np.polyfit(dff[3], hhh4, 1)
        slope55, intercept55 = np.polyfit(dff[4], hhh5, 1)
        slope66, intercept66 = np.polyfit(dff[5], hhh6, 1)

        yfit11 = [intercept11 + slope11 * xi for xi in dff[0]]
        yfit22 = [intercept22 + slope22 * xi for xi in dff[1]]
        yfit33 = [intercept33 + slope33 * xi for xi in dff[2]]
        yfit44 = [intercept44 + slope44 * xi for xi in dff[3]]
        yfit55 = [intercept55 + slope55 * xi for xi in dff[4]]
        yfit66 = [intercept66 + slope66 * xi for xi in dff[5]]

        slope=[slope1,slope2,slope3,slope4,slope5,slope6,0]
        intercept=[intercept1,intercept2,intercept3,intercept4,intercept5,intercept6,0]
        mx = [round(xi, 2) for xi in slope]
        
        slope1=[slope11,slope22,slope33,slope44,slope55,slope66,0]
        intercept1=[intercept11,intercept22,intercept33,intercept44,intercept55,intercept66,0]
        mx1 = [round(xi, 2) for xi in slope1]

        csvRow=[ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,dff,slope,intercept]
        csvfile = "/home/toyegoke/ENS_kinetics/Final-RDS-k/"+title+"_data_RDSk_kf_.csv"
        
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            #wr.writerow(csvRow)
            wr.writerow(['r_adsR','r_deh1','r_deh2','r_trnsV','r_desP','r_desH','ln(kf, P)'])
            wr.writerow([ppp1[0],ppp2[0],ppp3[0],ppp4[0],ppp5[0],ppp6[0],dff[0]])
            wr.writerow([ppp1[1],ppp2[1],ppp3[1],ppp4[1],ppp5[1],ppp6[1],dff[1]])
            wr.writerow([ppp1[2],ppp2[2],ppp3[2],ppp4[2],ppp5[2],ppp6[2],dff[2]])
            wr.writerow([ppp1[3],ppp2[3],ppp3[3],ppp4[3],ppp5[3],ppp6[3],dff[3]])
            wr.writerow([ppp1[4],ppp2[4],ppp3[4],ppp4[4],ppp5[4],ppp6[4],dff[4]])
            wr.writerow([ppp1[5],ppp2[5],ppp3[5],ppp4[5],ppp5[5],ppp6[5],dff[5]])
            wr.writerow([slope[0],slope[1],slope[2],slope[3],slope[4],slope[5],['slope']])
            wr.writerow([intercept[0],intercept[1],intercept[2],intercept[3],intercept[4],intercept[5],['intercept']])

        """
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(gaps,hhh1,marker='o',linestyle='--', label='adsR')
        ax.plot(gaps,hhh2,marker='x',linestyle='-.', label='deh1')
        ax.plot(gaps,hhh3,marker='v',linestyle='--', label='deh2')
        ax.plot(gaps,hhh4,marker='^',linestyle='-.', label='trns')
        ax.plot(gaps,hhh5,marker='+',linestyle='--', label='desP')
        ax.plot(gaps,hhh6,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G).png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(gaps,ppp1,marker='o',linestyle='--', label='adsR')
        ay.plot(gaps,ppp2,marker='x',linestyle='-.', label='deh1')
        ay.plot(gaps,ppp3,marker='v',linestyle='--', label='deh2')
        ay.plot(gaps,ppp4,marker='^',linestyle='-.', label='trns')
        ay.plot(gaps,ppp5,marker='+',linestyle='--', label='desP')
        ay.plot(gaps,ppp6,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G).png')
        """

        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(gaps,yfit11,marker='o',linestyle='--', label='adsR (m='+str(mx1[0])+')')
        ax.plot(gaps,yfit22,marker='x',linestyle='-.', label='deh1 (m='+str(mx1[1])+')')
        ax.plot(gaps,yfit33,marker='v',linestyle='--', label='deh2 (m='+str(mx1[2])+')')
        ax.plot(gaps,yfit44,marker='^',linestyle='-.', label='trns (m='+str(mx1[3])+')')
        ax.plot(gaps,yfit55,marker='+',linestyle='--', label='desP (m='+str(mx1[4])+')')
        ax.plot(gaps,yfit66,marker='x',linestyle='-.', label='desH (m='+str(mx1[5])+')')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G)_bestfitline.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(gaps,yfit1,marker='o',linestyle='--', label='adsR (m='+str(mx[0])+')')
        ax.plot(gaps,yfit2,marker='x',linestyle='-.', label='deh1 (m='+str(mx[1])+')')
        ax.plot(gaps,yfit3,marker='v',linestyle='--', label='deh2 (m='+str(mx[2])+')')
        ax.plot(gaps,yfit4,marker='^',linestyle='-.', label='trns (m='+str(mx[3])+')')
        ax.plot(gaps,yfit5,marker='+',linestyle='--', label='desP (m='+str(mx[4])+')')
        ax.plot(gaps,yfit6,marker='x',linestyle='-.', label='desH (m='+str(mx[5])+')')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G)_bestfitline.png')

        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(gaps,hhh1,marker='o',linestyle=' ')
        ax.plot(gaps,hhh2,marker='s',linestyle=' ')
        ax.plot(gaps,hhh3,marker='v',linestyle=' ')
        ax.plot(gaps,hhh4,marker='^',linestyle=' ')
        ax.plot(gaps,hhh5,marker='+',linestyle=' ')
        ax.plot(gaps,hhh6,marker='x',linestyle=' ')
        ax.plot(gaps,yfit11,linestyle='--', label='adsR (m='+str(mx1[0])+')')
        ax.plot(gaps,yfit22,linestyle='--', label='deh1 (m='+str(mx1[1])+')')
        ax.plot(gaps,yfit33,linestyle='--', label='deh2 (m='+str(mx1[2])+')')
        ax.plot(gaps,yfit44,linestyle='--', label='trns (m='+str(mx1[3])+')')
        ax.plot(gaps,yfit55,linestyle='--', label='desP (m='+str(mx1[4])+')')
        ax.plot(gaps,yfit66,linestyle='--', label='desH (m='+str(mx1[5])+')')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G)_scatterbestfitline.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(gaps,ppp1,marker='o',linestyle=' ')
        ay.plot(gaps,ppp2,marker='s',linestyle=' ')
        ay.plot(gaps,ppp3,marker='v',linestyle=' ')
        ay.plot(gaps,ppp4,marker='^',linestyle=' ')
        ay.plot(gaps,ppp5,marker='+',linestyle=' ')
        ay.plot(gaps,ppp6,marker='x',linestyle=' ')
        ay.plot(gaps,yfit1,linestyle='--', label='adsR (m='+str(mx[0])+')')
        ay.plot(gaps,yfit2,linestyle='--', label='deh1 (m='+str(mx[1])+')')
        ay.plot(gaps,yfit3,linestyle='--', label='deh2 (m='+str(mx[2])+')')
        ay.plot(gaps,yfit4,linestyle='--', label='trns (m='+str(mx[3])+')')
        ay.plot(gaps,yfit5,linestyle='--', label='desP (m='+str(mx[4])+')')
        ay.plot(gaps,yfit6,linestyle='--', label='desH (m='+str(mx[5])+')')
        plt.xlabel('$\Delta k_f$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G)_scatterbestfitline.png')

print('COMPLETED')
    
print('CALCULATION COMPLETED')

