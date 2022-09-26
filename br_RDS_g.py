from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdatab
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5,ppt6
from BR_kin_temp import kineticparameter, dydtb,solve_odes46b1
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
         #  BATCH - REACTOR
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def RDS(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout):      # <=============  

    #for specific time range GOOD FOR RDS

    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K   
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime ## in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp ## in K # temperature_selected # for each molecule study
    mc0 = cmass=1e-3 ## in mg
    P_R0 = 1.01325e5*10 ## Pa 
    C_R0 = P_R0/R/temp ## mol/m3 
    outflow_rate = flowout * (14.27) * 1e-7 ## in m3/sec
    v_R0 = (45*1e-6) ## in m3 (45 mL)
    n_R0 = v_R0*C_R0 ## in mol  

    ppt = moleculer_property_function # ppt1, ppt2

    S_x_in = (18*1e-6) ## in mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) ## in umol/m^2 (equivalent of 0.44 ML from expt data)
     
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    # Co = float(p_o/R/298.15) # in mol/m^3
    
    #Catalyst properties
    dc = 29 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2   
   
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-RDS-G/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-RDS-G/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)          


    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============    
        CCL = []; RR=[];  Rrr = []
        RP = []; RPP = []
        RH = []; RHH = []
        T = temp
            
        tt0, y= solve_odes46b1(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydtb, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)
 
 
        rxn_specie_id = Order_specieID_list     # list of selected specie IDs in a special order with string  elements
        x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
        vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; 
        ts2=rxn_specie_id[10]; ts3=rxn_specie_id[11]; tsr=rxn_specie_id[12]; tsp=rxn_specie_id[13]; tsh2=rxn_specie_id[14] 

        fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
        bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
        gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]; gat8=gat[7]; gat9=gat[8]; gat10=gat[9]; gat11=gat[10]
       
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
        TS3g=energylist[specie_index(ts3,int_id)]+gat8
        TSRg=energylist[specie_index(tsr,int_id)]+gat9
        TSPg=energylist[specie_index(tsp,int_id)]+gat10
        TSH2g=energylist[specie_index(tsh2,int_id)]+gat11 
       
        # calculate all reaction rate constants
        print('2R+X<==>RX')
        kf1=kineticparameter.ksu(2*Rg+Xg,2*Rg+Xg,S_o,2,1,T)*fat1 
        kb1=kineticparameter.ksu(2*Rg+Xg,RXg,S_o,0,1,T)*bat1 
        """kf1=kineticparameter.ksu(TSRg,2*Rg+Xg,S_o,2,1,T)*fat1 
        kb1=kineticparameter.ksu(TSRg,RXg,S_o,0,1,T)*bat1 """

        print('RX+X<==>UX+HX')
        kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)*fat2 
        kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)*bat2

        print('UX+X<==>VX+HX')
        kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)*fat3 
        kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)*bat3

        print('VX+X<==>PX+HX')
        kf6=kineticparameter.ksu(TS3g+Xg,VXg+Xg,S_o,0,2,T)*fat4 
        kb6=kineticparameter.ksu(TS3g+Xg,PXg+HXg,S_o,0,2,T)*bat4  

        print('PX<==>P+X')
        kb4=kineticparameter.ksu(Pg+Xg,Pg+Xg,S_o,1,1,T)*bat5  
        kf4=kineticparameter.ksu(Pg+Xg,PXg,S_o,0,1,T)*fat5  
        """kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T)*bat5  
        kf4=kineticparameter.ksu(TSPg,PXg,S_o,0,1,T)*fat5"""  
    
        print('H2X<==>H2+X')
        kb5=kineticparameter.ksu(H2g+Xg,H2g+Xg,S_o,1,1,T)*bat6 
        kf5=kineticparameter.ksu(H2g+Xg,HXg,S_o,0,1,T)*fat6 
        """kb5=kineticparameter.ksu(TSH2g,H2g+Xg,S_o,1,1,T)*bat6 
        kf5=kineticparameter.ksu(TSH2g,HXg,S_o,0,1,T)*fat6""" 
 
        
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

        # RR.append(rr)

        RP.append(rP)
        RH.append(rH)

        gatG = [RXg,UXg,HXg,VXg,PXg,TS1g,TS2g,TSRg,TSPg,TSH2g]
        print(gatG)    

        Rrr.append((n_R/n_Ro - n_Rb/n_Ro)/(Ttim - Ttimb)) 
        RPP.append((n_P/n_Ro - n_Pb/n_Ro)/(Ttim - Ttimb)) 
        RHH.append((n_H2/n_Ro - n_H2b/n_Ro)/(Ttim - Ttimb))

        # Rrr.append(-outflow_rate/v_R0*(n_Ro-n_R)/n_Ro) # in 1/s 
        # RPP.append(-outflow_rate/v_R0*(0-n_P)/n_Ro) # in 1/s
        # RHH.append(-outflow_rate/v_R0*(0-n_H2)/n_Ro)  # in 1/s

        # return RH[0], RP[0], loge(RH[0]), loge(RP[0]), gatG    
        return RHH[0], RPP[0], loge(RHH[0]), loge(RPP[0]), gatG     
  

    def gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list):
    
        # title0=titlemolecule
        dff = [-91000,-51000,-21000,0,21000,51000,91000]

        #gat = [0,0,0,0,0,0,0]
        bat = [1,1,1,1,1,1]
        fat = bat 
        T = temp

        ppt = moleculer_property_function # ppt1, ppt2,

        rangedff=len(dff)
        
        rrr1=[]; ppp1=[]; hhh1=[]; rrr2=[]; ppp2=[]; hhh2=[] 
        rrr3=[]; ppp3=[]; hhh3=[]; rrr4=[]; ppp4=[]; hhh4=[] 
        rrr5=[]; ppp5=[]; hhh5=[]; rrr6=[]; ppp6=[]; hhh6=[]
        rrr7=[]; ppp7=[]; hhh7=[]; rrr8=[]; ppp8=[]; hhh8=[]
        rrr9=[]; ppp9=[]; hhh9=[]; rrr10=[]; ppp10=[]; hhh10=[]
        rrr11=[]; ppp11=[]; hhh11=[]
     
        i=0    
        lnG1=[]; lnG2=[]; lnG3=[]; lnG4=[]; lnG5=[]; lnG6=[]; lnG7=[]; 
        lnG8=[]; lnG9=[]; lnG10=[]; lnG11=[]

        for i in range(rangedff):
        
            gat1=[dff[i],0,0,0,0,0,0,0, 0,0,0] 
            gat2=[0,dff[i],0,0,0,0,0,0, 0,0,0] 
            gat3=[0,0,dff[i],0,0,0,0,0, 0,0,0] 
            gat4=[0,0,0,dff[i],0,0,0,0, 0,0,0] 
            gat5=[0,0,0,0,dff[i],0,0,0, 0,0,0] 
            gat6=[0,0,0,0,0,dff[i],0,0, 0,0,0] 
            gat7=[0,0,0,0,0,0,dff[i],0, 0,0,0] 
            gat8=[0,0,0,0,0,0,0,dff[i], 0,0,0] 
            gat9=[0,0,0,0,0,0,0,0,dff[i], 0,0] 
            gat10=[0,0,0,0,0,0,0,0,0,dff[i],0] 
            gat11=[0,0,0,0,0,0,0,0,0,0,dff[i]] 
        
            # the arguments are ====>   mc0,  C_R0,  outflow_rate,  Temp, time, dydt,  Gm, intID, specieID_list, S_x, fat, bat, gat
            RH1, RP1, lnRH1, lnRP1, gatG1 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat1)
            RH2, RP2, lnRH2, lnRP2, gatG2 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat2)
            RH3, RP3, lnRH3, lnRP3, gatG3 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat3)
            RH4, RP4, lnRH4, lnRP4, gatG4 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat4)
            RH5, RP5, lnRH5, lnRP5, gatG5 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat5)
            RH6, RP6, lnRH6, lnRP6, gatG6 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat6)
            RH7, RP7, lnRH7, lnRP7, gatG7 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat7)
            RH8, RP8, lnRH8, lnRP8, gatG8 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat8)
            RH9, RP9, lnRH9, lnRP9, gatG9 = RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat9)
            RH10,RP10,lnRH10,lnRP10,gatG10= RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat10)
            RH11,RP11,lnRH11,lnRP11,gatG11= RDS_effect2(ppt, title0, Order_specieID_list, fat, bat, gat11)

            hhh1.append(lnRH1); ppp1.append(lnRP1) 
            hhh2.append(lnRH2); ppp2.append(lnRP2) 
            hhh3.append(lnRH3); ppp3.append(lnRP3) 
            hhh4.append(lnRH4); ppp4.append(lnRP4) 
            hhh5.append(lnRH5); ppp5.append(lnRP5) 
            hhh6.append(lnRH6); ppp6.append(lnRP6)
            hhh7.append(lnRH7); ppp7.append(lnRP7)
            hhh8.append(lnRH8); ppp8.append(lnRP8)
            hhh9.append(lnRH9); ppp9.append(lnRP9)
            hhh10.append(lnRH10); ppp10.append(lnRP10)
            hhh11.append(lnRH11); ppp11.append(lnRP11)

            lnG1.append(float(-gatG1[0]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG2.append(float(-gatG2[1]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG3.append(float(-gatG3[2]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG4.append(float(-gatG4[3]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG5.append(float(-gatG5[4]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG6.append(float(-gatG6[5]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG7.append(float(-gatG7[6]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG8.append(float(-gatG8[7]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG9.append(float(-gatG9[8]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG10.append(float(-gatG10[9]/R/T)) # the nat log of the change in G i.e. ln(G)
            lnG11.append(float(-gatG11[-1]/R/T)) # the nat log of the change in G i.e. ln(G)

        lnG=[lnG1,lnG2,lnG3,lnG4,lnG5,lnG6,lnG7,lnG8,lnG9,lnG10,lnG11]; gaps=[dff,dff,dff,dff,dff,dff,dff,dff,dff,dff]
        print(lnG)

        return hhh1,hhh2,hhh3,hhh4,hhh5,hhh6,hhh7,hhh8,hhh9,hhh10,hhh11,ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,ppp7,ppp8,ppp9,ppp10,ppp11,lnG,gaps

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

title1= 'Ammonia Dehydrogenation over Ru'
engg1 = ['X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TS3','TSR','TSP','TSH2']
TIM=3.39
pt = [ppt6]
titl = [title1]
specie_ID = [engg1]
cmass = 50 # 1e-3 # mg
flowout = 0 # semi-batch is 1 and batch is 0
rtemp = [700+273.15] # in K
rtime = [10**(5)] # in sec
ranP = len(titl)
i=0

for i in range(ranP):
        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        hhh1,hhh2,hhh3,hhh4,hhh5,hhh6,hhh7,hhh8,hhh9,hhh10,hhh11,ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,ppp7,ppp8,ppp9,ppp10,ppp11,dff,gaps=RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout)
        
        slope1, intercept1 = np.polyfit(dff[0], ppp1, 1)
        slope2, intercept2 = np.polyfit(dff[1], ppp2, 1)
        slope3, intercept3 = np.polyfit(dff[2], ppp3, 1)
        slope4, intercept4 = np.polyfit(dff[3], ppp4, 1)
        slope5, intercept5 = np.polyfit(dff[4], ppp5, 1)
        slope6, intercept6 = np.polyfit(dff[5], ppp6, 1)
        slope7, intercept7 = np.polyfit(dff[6], ppp7, 1)
        slope8, intercept8 = np.polyfit(dff[7], ppp8, 1)
        slope9, intercept9 = np.polyfit(dff[8], ppp9, 1)
        slope10, intercept10 = np.polyfit(dff[9], ppp10, 1)
        slope11, intercept11 = np.polyfit(dff[10], ppp11, 1)

        yfit1 = [intercept1 + slope1 * float(xi) for xi in dff[0]]
        yfit2 = [intercept2 + slope2 * float(xi) for xi in dff[1]]
        yfit3 = [intercept3 + slope3 * float(xi) for xi in dff[2]]
        yfit4 = [intercept4 + slope4 * float(xi) for xi in dff[3]]
        yfit5 = [intercept5 + slope5 * float(xi) for xi in dff[4]]
        yfit6 = [intercept6 + slope6 * float(xi) for xi in dff[5]]
        yfit7 = [intercept7 + slope7 * float(xi) for xi in dff[6]]
        yfit8 = [intercept8 + slope8 * float(xi) for xi in dff[7]]
        yfit9 = [intercept9 + slope9 * float(xi) for xi in dff[8]]
        yfit10 = [intercept10 + slope10 * float(xi) for xi in dff[9]]
        yfit11 = [intercept11 + slope11 * float(xi) for xi in dff[10]]

        slope11, intercept11 = np.polyfit(dff[0], hhh1, 1)
        slope22, intercept22 = np.polyfit(dff[1], hhh2, 1)
        slope33, intercept33 = np.polyfit(dff[2], hhh3, 1)
        slope44, intercept44 = np.polyfit(dff[3], hhh4, 1)
        slope55, intercept55 = np.polyfit(dff[4], hhh5, 1)
        slope66, intercept66 = np.polyfit(dff[5], hhh6, 1)
        slope77, intercept77 = np.polyfit(dff[6], hhh7, 1)
        slope88, intercept88 = np.polyfit(dff[7], hhh8, 1)
        slope99, intercept99 = np.polyfit(dff[8], hhh9, 1)
        slope100, intercept100 = np.polyfit(dff[9], hhh10, 1)
        slope110, intercept110 = np.polyfit(dff[10], hhh11, 1)

        yfit11 = [intercept11 + slope11 * xi for xi in dff[0]]
        yfit22 = [intercept22 + slope22 * xi for xi in dff[1]]
        yfit33 = [intercept33 + slope33 * xi for xi in dff[2]]
        yfit44 = [intercept44 + slope44 * xi for xi in dff[3]]
        yfit55 = [intercept55 + slope55 * xi for xi in dff[4]]
        yfit66 = [intercept66 + slope66 * xi for xi in dff[5]]
        yfit77 = [intercept77 + slope77 * xi for xi in dff[6]]
        yfit88 = [intercept88 + slope88 * xi for xi in dff[7]]
        yfit99 = [intercept99 + slope99 * xi for xi in dff[8]]
        yfit100 = [intercept100 + slope100 * xi for xi in dff[9]]
        yfit110 = [intercept110 + slope110 * xi for xi in dff[10]]

        slope=[slope1,slope2,slope3,slope4,slope5,slope6,slope7,slope8,slope9,slope10,slope11,0]
        intercept=[intercept1,intercept2,intercept3,intercept4,intercept5,intercept6,intercept7,intercept8,intercept9,intercept10,intercept11,0]
        mx = [round(xi, 2) for xi in slope]
        
        slope1=[slope11,slope22,slope33,slope44,slope55,slope66,slope77,slope88,slope99,slope100,slope110,0]
        intercept1=[intercept11,intercept22,intercept33,intercept44,intercept55,intercept66,intercept77,intercept88,intercept99,intercept100,intercept110,0]
        mx1 = [round(xi, 2) for xi in slope1]

        csvRow=[ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,ppp7,ppp8,ppp9,ppp10,ppp11,dff,slope,intercept]
        csvfile = "/home/toyegoke/ENS_kinetics/Final-RDS-G/"+title+"_data_RDS_G_.csv"
        
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            #wr.writerow(csvRow)
            wr.writerow(['r_RX','r_UX','r_HX','r_VX','r_PX','r_TS1','r_TS2','ln(G,P)'])
            wr.writerow([ppp1[0],ppp2[0],ppp3[0],ppp4[0],ppp5[0],ppp6[0],ppp7[0],dff[0]])
            wr.writerow([ppp1[1],ppp2[1],ppp3[1],ppp4[1],ppp5[1],ppp6[1],ppp7[1],dff[1]])
            wr.writerow([ppp1[2],ppp2[2],ppp3[2],ppp4[2],ppp5[2],ppp6[2],ppp7[2],dff[2]])
            wr.writerow([ppp1[3],ppp2[3],ppp3[3],ppp4[3],ppp5[3],ppp6[3],ppp7[3],dff[3]])
            wr.writerow([ppp1[4],ppp2[4],ppp3[4],ppp4[4],ppp5[4],ppp6[4],ppp7[4],dff[4]])
            wr.writerow([ppp1[5],ppp2[5],ppp3[5],ppp4[5],ppp5[5],ppp6[5],ppp7[5],dff[5]])
            wr.writerow([ppp1[6],ppp2[6],ppp3[6],ppp4[6],ppp5[6],ppp6[6],ppp7[6],dff[6]])
            wr.writerow([slope[0],slope[1],slope[2],slope[3],slope[4],slope[5],slope[6],['slope']])
            wr.writerow([intercept[0],intercept[1],intercept[2],intercept[3],intercept[4],intercept[5],intercept[6],['intercept']])

        """
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(hhh1,marker='o',linestyle='--', label='RX')
        ax.plot(dff[1],hhh2,marker='x',linestyle='-.', label='UX')
        ax.plot(dff[2],hhh3,marker='v',linestyle='--', label='HX')
        ax.plot(dff[3],hhh4,marker='^',linestyle='-.', label='VX')
        ax.plot(dff[4],hhh5,marker='+',linestyle='--', label='PX')
        ax.plot(dff[5],hhh6,marker='x',linestyle='-.', label='TS1')
        ax.plot(dff[6],hhh7,marker='*',linestyle='-.', label='TS2')
        ax.plot(dff[7],hhh8,marker='*',linestyle='-.', label='TSR')
        ax.plot(dff[8],hhh9,marker='*',linestyle='-.', label='TSP')
        ax.plot(dff[9],hhh0,marker='*',linestyle='-.', label='TSH2')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(H2_G).png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(dff[0],ppp1,marker='o',linestyle='--', label='RX')
        ay.plot(dff[1],ppp2,marker='x',linestyle='-.', label='UX')
        ay.plot(dff[2],ppp3,marker='v',linestyle='--', label='HX')
        ay.plot(dff[3],ppp4,marker='^',linestyle='-.', label='VX')
        ay.plot(dff[4],ppp5,marker='+',linestyle='--', label='PX')
        ay.plot(dff[5],ppp6,marker='x',linestyle='-.', label='TS1')
        ay.plot(dff[6],ppp7,marker='*',linestyle='-.', label='TS2')
        ay.plot(dff[7],ppp8,marker='*',linestyle='-.', label='TSR')
        ay.plot(dff[8],ppp9,marker='*',linestyle='-.', label='TSP')
        ay.plot(dff[9],ppp0,marker='*',linestyle='-.', label='TSH2')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(P_G).png')
        """


        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],yfit11,marker='o',linestyle='--', label='RX (m='+str(mx1[0])+')')
        ax.plot(dff[1],yfit22,marker='x',linestyle='-.', label='UX (m='+str(mx1[1])+')')
        ax.plot(dff[2],yfit33,marker='v',linestyle='--', label='HX (m='+str(mx1[2])+')')
        ax.plot(dff[3],yfit44,marker='^',linestyle='-.', label='VX (m='+str(mx1[3])+')')
        ax.plot(dff[4],yfit55,marker='+',linestyle='--', label='PX (m='+str(mx1[4])+')')
        ax.plot(dff[5],yfit66,marker='x',linestyle='-.', label='TS1 (m='+str(mx1[5])+')')
        ax.plot(dff[6],yfit77,marker='x',linestyle='-.', label='TS2 (m='+str(mx1[6])+')')
        ax.plot(dff[7],yfit88,marker='x',linestyle='-.', label='TS3 (m='+str(mx1[7])+')')
        ax.plot(dff[8],yfit99,marker='x',linestyle='-.', label='TSR (m='+str(mx1[8])+')')
        ax.plot(dff[9],yfit100,marker='x',linestyle='-.', label='TSP (m='+str(mx1[9])+')')
        ax.plot(dff[10],yfit110,marker='x',linestyle='-.', label='TSH2 (m='+str(mx1[10])+')')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(H2_G)_bestfitline.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],yfit1,marker='o',linestyle='--', label='RX (m='+str(mx[0])+')')
        ax.plot(dff[1],yfit2,marker='x',linestyle='-.', label='UX (m='+str(mx[1])+')')
        ax.plot(dff[2],yfit3,marker='v',linestyle='--', label='HX (m='+str(mx[2])+')')
        ax.plot(dff[3],yfit4,marker='^',linestyle='-.', label='VX (m='+str(mx[3])+')')
        ax.plot(dff[4],yfit5,marker='+',linestyle='--', label='PX (m='+str(mx[4])+')')
        ax.plot(dff[5],yfit6,marker='x',linestyle='-.', label='TS1 (m='+str(mx[5])+')')
        ax.plot(dff[6],yfit7,marker='x',linestyle='-.', label='TS2 (m='+str(mx[6])+')')
        ax.plot(dff[7],yfit8,marker='x',linestyle='-.', label='TS3 (m='+str(mx[7])+')')
        ax.plot(dff[8],yfit9,marker='x',linestyle='-.', label='TSR (m='+str(mx[8])+')')
        ax.plot(dff[9],yfit10,marker='x',linestyle='-.', label='TSP (m='+str(mx[9])+')')
        ax.plot(dff[10],yfit11,marker='x',linestyle='-.', label='TSH2 (m='+str(mx[10])+')')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(P_G)_bestfitline.png')

        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.scatter(dff[0],hhh1,marker='o')
        ax.scatter(dff[1],hhh2,marker='s')
        ax.scatter(dff[2],hhh3,marker='v')
        ax.scatter(dff[3],hhh4,marker='^')
        ax.scatter(dff[4],hhh5,marker='+')
        ax.scatter(dff[5],hhh6,marker='x')
        ax.scatter(dff[6],hhh7,marker='*')
        ax.scatter(dff[7],hhh8,marker='*')
        ax.scatter(dff[8],hhh9,marker='*')
        ax.scatter(dff[9],hhh10,marker='*')
        ax.scatter(dff[10],hhh11,marker='*')
        ax.plot(dff[0],yfit11,linestyle='--', label='RX (m='+str(mx1[0])+')')
        ax.plot(dff[1],yfit22,linestyle='-.', label='UX (m='+str(mx1[1])+')')
        ax.plot(dff[2],yfit33,linestyle='--', label='HX (m='+str(mx1[2])+')')
        ax.plot(dff[3],yfit44,linestyle='-.', label='VX (m='+str(mx1[3])+')')
        ax.plot(dff[4],yfit55,linestyle='--', label='PX (m='+str(mx1[4])+')')
        ax.plot(dff[5],yfit66,linestyle='-.', label='TS1 (m='+str(mx1[5])+')')
        ax.plot(dff[6],yfit77,linestyle='-.', label='TS2 (m='+str(mx1[6])+')')
        ax.plot(dff[7],yfit88,linestyle='-.', label='TS3 (m='+str(mx1[7])+')')
        ax.plot(dff[8],yfit99,linestyle='-.', label='TSR (m='+str(mx1[8])+')')
        ax.plot(dff[9],yfit100,linestyle='-.', label='TSP (m='+str(mx1[9])+')')
        ax.plot(dff[10],yfit110,linestyle='-.', label='TSH2 (m='+str(mx1[10])+')')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(H2_G)_scatterbestfitline.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.scatter(dff[0],ppp1,marker='o')
        ay.scatter(dff[1],ppp2,marker='s')
        ay.scatter(dff[2],ppp3,marker='v')
        ay.scatter(dff[3],ppp4,marker='^')
        ay.scatter(dff[4],ppp5,marker='+')
        ay.scatter(dff[5],ppp6,marker='x')
        ay.scatter(dff[6],ppp7,marker='*')
        ay.scatter(dff[7],ppp8,marker='*')
        ay.scatter(dff[8],ppp9,marker='*')
        ay.scatter(dff[9],ppp10,marker='*')
        ay.scatter(dff[10],ppp11,marker='*')
        ay.plot(dff[0],yfit1,linestyle='--', label='RX (m='+str(mx[0])+')')
        ay.plot(dff[1],yfit2,linestyle='-.', label='UX (m='+str(mx[1])+')')
        ay.plot(dff[2],yfit3,linestyle='--', label='HX (m='+str(mx[2])+')')
        ay.plot(dff[3],yfit4,linestyle='-.', label='VX (m='+str(mx[3])+')')
        ay.plot(dff[4],yfit5,linestyle='--', label='PX (m='+str(mx[4])+')')
        ay.plot(dff[5],yfit6,linestyle='-.', label='TS1(m='+str(mx[5])+')')
        ay.plot(dff[6],yfit7,linestyle='-.', label='TS2(m='+str(mx[6])+')')
        ay.plot(dff[7],yfit8,linestyle='-.', label='TS3(m='+str(mx[7])+')')
        ay.plot(dff[8],yfit9,linestyle='-.', label='TSR(m='+str(mx[8])+')')
        ay.plot(dff[9],yfit10,linestyle='-.', label='TSP(m='+str(mx[9])+')')
        ay.plot(dff[10],yfit11,linestyle='-.', label='TSH2(m='+str(mx[10])+')')
        plt.xlabel('$(^{-G}/_{RT})$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-G/'+title+'_plot_data_for_RDS(P_G)_scatterbestfitline.png')

print(ppp1, dff)
print('COMPLETED')
    
print('CALCULATION COMPLETED')






