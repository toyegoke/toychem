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

#===================================================================================================================================    

# TYPE 2 WITH TEMPERATURE EFFECT    

#===================================================================================================================================    

def RDStemp(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout):      # <=============  

    #for specific time range GOOD FOR RDS
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime # (14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp # 550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = cmass # 40  # in g
    C_R0 = 997.74 # mol/m3 
    outflow_rate =  flowout * (14.27) * 1e-7 # in m3/sec
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
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)          

    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============  
  
        CCL = []        #RR = []
        RP = []
        RH = []
        T = temp
            
        tt0, y= solve_odes40(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)

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
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_'+str(int(T))+'K_'+'_1amount_plot_.png')

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
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_'+str(int(T))+'K_'+'_1surface_plot_.png')
        #========================================

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

        RP.append(rP)
        RH.append(rH)

        fatK = [kf1,kf2,kf3,kf6,kf4,kf5]
        print(fatK)    

        return RH[0], RP[0], loge(RH[0]), loge(RP[0]), fatK       

    def gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list, temp):
    
        # title0=titlemolecule
        dff = [0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9]

        gat = [0,0,0,0,0,0,0]           #bat = [1,1,1,1,1,1]         #fat = bat 
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

    def RDStempsolve(moleculer_property_function,title_of_molecule,Order_specieID_list,temp):
        TT = [temp-20,temp-16,temp-12,temp-8,temp-4,temp,temp+4,temp+8,temp+12,temp+16,temp+20]
        m1=[]; m2=[]; m3=[]; m4=[]; m5=[]; m6=[]
        m11=[]; m22=[]; m33=[]; m44=[]; m55=[]; m66=[]
        c1=[]; c2=[]; c3=[]; c4=[]; c5=[]; c6=[]
        c11=[]; c22=[]; c33=[]; c44=[]; c55=[]; c66=[]  
        tlens = len(TT)
        i = 0
        for i in range(tlens):
            T=TT[i]
            hhh1,hhh2,hhh3,hhh4,hhh5,hhh6,ppp1,ppp2,ppp3,ppp4,ppp5,ppp6,dff,gaps =gat_RDS2(ppt,title,engg,T)
            
            slope1, intercept1 = np.polyfit(dff[0], ppp1, 1)
            slope2, intercept2 = np.polyfit(dff[1], ppp2, 1)
            slope3, intercept3 = np.polyfit(dff[2], ppp3, 1)
            slope4, intercept4 = np.polyfit(dff[3], ppp4, 1)
            slope5, intercept5 = np.polyfit(dff[4], ppp5, 1)
            slope6, intercept6 = np.polyfit(dff[5], ppp6, 1)
    
            slope11, intercept11 = np.polyfit(dff[0], hhh1, 1)
            slope22, intercept22 = np.polyfit(dff[1], hhh2, 1)
            slope33, intercept33 = np.polyfit(dff[2], hhh3, 1)
            slope44, intercept44 = np.polyfit(dff[3], hhh4, 1)
            slope55, intercept55 = np.polyfit(dff[4], hhh5, 1)
            slope66, intercept66 = np.polyfit(dff[5], hhh6, 1)
    
            m1.append(slope1); m2.append(slope2); m3.append(slope3); m4.append(slope4); m5.append(slope5); m6.append(slope6)
            c1.append(intercept1); c2.append(intercept2); c3.append(intercept3); c4.append(intercept4); c5.append(intercept5); c6.append(intercept6)
            
            m11.append(slope11); m22.append(slope22); m33.append(slope33); m44.append(slope44); m55.append(slope55); m66.append(slope66)
            c11.append(intercept11); c22.append(intercept22); c33.append(intercept33); c44.append(intercept44)
            c55.append(intercept55); c66.append(intercept66)
    
        fig = plt.figure()
        ax = plt.subplot(111)
        gaps = TT
        ax.plot(gaps,m11,marker='o',linestyle='--', label='adsR')
        ax.plot(gaps,m22,marker='x',linestyle='-.', label='deh1')
        ax.plot(gaps,m33,marker='v',linestyle='--', label='deh2')
        ax.plot(gaps,m44,marker='^',linestyle='-.', label='trns')
        ax.plot(gaps,m55,marker='+',linestyle='--', label='desP')
        ax.plot(gaps,m66,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$T / K$')
        plt.ylabel('$X_{DRC}$')
        plt.title('DRC for '+title+' Kinetics')
        ax.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_plot_data_for_DRC_ktemp(H2_G).png')
    
        fig = plt.figure()
        ay = plt.subplot(111)
        gaps = TT
        ay.plot(gaps,m1,marker='o',linestyle='--', label='adsR')
        ay.plot(gaps,m2,marker='x',linestyle='-.', label='deh1')
        ay.plot(gaps,m3,marker='v',linestyle='--', label='deh2')
        ay.plot(gaps,m4,marker='^',linestyle='-.', label='trns')
        ay.plot(gaps,m5,marker='+',linestyle='--', label='desP')
        ay.plot(gaps,m6,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$T / K$')
        plt.ylabel('$X_{DRC}$')
        plt.title('DRC for '+title+' Kinetics')
        ay.legend()
        #plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/'+title+'_plot_data_for_DRC_ktemp(P_G).png')
    
        csvfile = "/home/toyegoke/ENS_kinetics/Final-RDS-ktemp/"+title+"_data_DRC_kf_TEMP_.csv"
            
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerow(['X_adsR','X_deh1','X_deh2','X_trnsV','X_desP','X_desH','mH2_T / K)'])
            wr.writerow([m11[0],m22[0],m33[0],m44[0],m55[0],m66[0],TT[0]])
            wr.writerow([m11[1],m22[1],m33[1],m44[1],m55[1],m66[1],TT[1]])
            wr.writerow([m11[2],m22[2],m33[2],m44[2],m55[2],m66[2],TT[2]])
            wr.writerow([m11[3],m22[3],m33[3],m44[3],m55[3],m66[3],TT[3]])
            wr.writerow([m11[4],m22[4],m33[4],m44[4],m55[4],m66[4],TT[4]])
            wr.writerow([m11[5],m22[5],m33[5],m44[5],m55[5],m66[5],TT[5]])
            wr.writerow([m11[6],m22[6],m33[6],m44[6],m55[6],m66[6],TT[6]])
            wr.writerow([m11[7],m22[7],m33[7],m44[7],m55[7],m66[7],TT[7]])
            wr.writerow([m11[8],m22[8],m33[8],m44[8],m55[8],m66[8],TT[8]])
            wr.writerow([m11[9],m22[9],m33[9],m44[9],m55[9],m66[9],TT[9]])
            wr.writerow([m11[10],m22[10],m3[10],m4[10],m5[10],m6[10],TT[10]])
            wr.writerow(['X_adsR','X_deh1','X_deh2','X_trnsV','X_desP','X_desH','cH2_T / K)'])
            wr.writerow([c11[0],c22[0],c33[0],c44[0],c55[0],c66[0],TT[0]])
            wr.writerow([c11[1],c22[1],c33[1],c44[1],c55[1],c66[1],TT[1]])
            wr.writerow([c11[2],c22[2],c33[2],c44[2],c55[2],c66[2],TT[2]])
            wr.writerow([c11[3],c22[3],c33[3],c44[3],c55[3],c66[3],TT[3]])
            wr.writerow([c11[4],c22[4],c33[4],c44[4],c55[4],c66[4],TT[4]])
            wr.writerow([c11[5],c22[5],c33[5],c44[5],c55[5],c66[5],TT[5]])
            wr.writerow([c11[6],c22[6],c33[6],c44[6],c55[6],c66[6],TT[6]])
            wr.writerow([c11[7],c22[7],c33[7],c44[7],c55[7],c66[7],TT[7]])
            wr.writerow([c11[8],c22[8],c33[8],c44[8],c55[8],c66[8],TT[8]])
            wr.writerow([c11[9],c22[9],c33[9],c44[9],c55[9],c66[9],TT[9]])
            wr.writerow([c11[10],c22[10],c33[10],c44[10],c55[10],c66[10],TT[10]])
    
            wr.writerow(['X_adsR','X_deh1','X_deh2','X_trnsV','X_desP','X_desH','mP_T / K)'])
            wr.writerow([m1[0],m2[0],m3[0],m4[0],m5[0],m6[0],TT[0]])
            wr.writerow([m1[1],m2[1],m3[1],m4[1],m5[1],m6[1],TT[1]])
            wr.writerow([m1[2],m2[2],m3[2],m4[2],m5[2],m6[2],TT[2]])
            wr.writerow([m1[3],m2[3],m3[3],m4[3],m5[3],m6[3],TT[3]])
            wr.writerow([m1[4],m2[4],m3[4],m4[4],m5[4],m6[4],TT[4]])
            wr.writerow([m1[5],m2[5],m3[5],m4[5],m5[5],m6[5],TT[5]])
            wr.writerow([m1[6],m2[6],m3[6],m4[6],m5[6],m6[6],TT[6]])
            wr.writerow([m1[7],m2[7],m3[7],m4[7],m5[7],m6[7],TT[7]])
            wr.writerow([m1[8],m2[8],m3[8],m4[8],m5[8],m6[8],TT[8]])
            wr.writerow([m1[9],m2[9],m3[9],m4[9],m5[9],m6[9],TT[9]])
            wr.writerow([m1[10],m2[10],m3[10],m4[10],m5[10],m6[10],TT[10]])
            wr.writerow(['X_adsR','X_deh1','X_deh2','X_trnsV','X_desP','X_desH','cP_T / K)'])
            wr.writerow([c1[0],c2[0],c3[0],c4[0],c5[0],c6[0],TT[0]])
            wr.writerow([c1[1],c2[1],c3[1],c4[1],c5[1],c6[1],TT[1]])
            wr.writerow([c1[2],c2[2],c3[2],c4[2],c5[2],c6[2],TT[2]])
            wr.writerow([c1[3],c2[3],c3[3],c4[3],c5[3],c6[3],TT[3]])
            wr.writerow([c1[4],c2[4],c3[4],c4[4],c5[4],c6[4],TT[4]])
            wr.writerow([c1[5],c2[5],c3[5],c4[5],c5[5],c6[5],TT[5]])
            wr.writerow([c1[6],c2[6],c3[6],c4[6],c5[6],c6[6],TT[6]])
            wr.writerow([c1[7],c2[7],c3[7],c4[7],c5[7],c6[7],TT[7]])
            wr.writerow([c1[8],c2[8],c3[8],c4[8],c5[8],c6[8],TT[8]])
            wr.writerow([c1[9],c2[9],c3[9],c4[9],c5[9],c6[9],TT[9]])
            wr.writerow([c1[10],c2[10],c3[10],c4[10],c5[10],c6[10],TT[10]])
    
        return csvfile

    return RDStempsolve(moleculer_property_function,title_of_molecule,Order_specieID_list,temp)


#=============================================================================================================================


    
def best_fit(X, Y):
    # solve for a and b // intercept & slope
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))
    
    return a, b


#===================================================================================================================================

# call for solution section

#===================================================================================================================================


#================================================
# for RDS with temperature effect
#================================================

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

rtime = [3854.608589,666.2822896,3214.544089,2375.114529,1377.532526] # for Tm+25K +5%
#rtime = [30181.24358,36190.7871,5843905.391,18596.93624,2088453.228] # for Tm+25K +50%
rtemp = [502,460,418,502,539] # Tmin+25K
#rtime = [59.18063593,43.72650814,55.70480542,32.30799204,555.6449496] # for 550K +5%
#rtime = [463.3791336,798.9492025,901.7644281,303.3386229,702521.0872] # for 550K +50%
#rtemp = [550,550,550,550,550] # 550K

cmass = 1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

ranP = len(titl)
i=0

for i in range(ranP):

        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        solution =RDStemp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout)

print('COMPLETED')
    
print('CALCULATION COMPLETED')






































#3333333333333333333333333333333333333333333333333333333

"""
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],hhh1,marker='o',linestyle='--', label='adsR')
        ax.plot(dff[1],hhh2,marker='x',linestyle='-.', label='deh1')
        ax.plot(dff[2],hhh3,marker='v',linestyle='--', label='deh2')
        ax.plot(dff[3],hhh4,marker='^',linestyle='-.', label='trns')
        ax.plot(dff[4],hhh5,marker='+',linestyle='--', label='desP')
        ax.plot(dff[5],hhh6,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G).png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(dff[0],ppp1,marker='o',linestyle='--', label='adsR')
        ay.plot(dff[1],ppp2,marker='x',linestyle='-.', label='deh1')
        ay.plot(dff[2],ppp3,marker='v',linestyle='--', label='deh2')
        ay.plot(dff[3],ppp4,marker='^',linestyle='-.', label='trns')
        ay.plot(dff[4],ppp5,marker='+',linestyle='--', label='desP')
        ay.plot(dff[5],ppp6,marker='x',linestyle='-.', label='desH')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G).png')

        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],yfit11,marker='o',linestyle='--', label='adsR (m='+str(mx1[0])+')')
        ax.plot(dff[1],yfit22,marker='x',linestyle='-.', label='deh1 (m='+str(mx1[1])+')')
        ax.plot(dff[2],yfit33,marker='v',linestyle='--', label='deh2 (m='+str(mx1[2])+')')
        ax.plot(dff[3],yfit44,marker='^',linestyle='-.', label='trns (m='+str(mx1[3])+')')
        ax.plot(dff[4],yfit55,marker='+',linestyle='--', label='desP (m='+str(mx1[4])+')')
        ax.plot(dff[5],yfit66,marker='x',linestyle='-.', label='desH (m='+str(mx1[5])+')')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G)_bestfitline.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],yfit1,marker='o',linestyle='--', label='adsR (m='+str(mx[0])+')')
        ax.plot(dff[1],yfit2,marker='x',linestyle='-.', label='deh1 (m='+str(mx[1])+')')
        ax.plot(dff[2],yfit3,marker='v',linestyle='--', label='deh2 (m='+str(mx[2])+')')
        ax.plot(dff[3],yfit4,marker='^',linestyle='-.', label='trns (m='+str(mx[3])+')')
        ax.plot(dff[4],yfit5,marker='+',linestyle='--', label='desP (m='+str(mx[4])+')')
        ax.plot(dff[5],yfit6,marker='x',linestyle='-.', label='desH (m='+str(mx[5])+')')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G)_bestfitline.png')

        #xxxxxxxxxxxxxxxxxxxxxx

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(dff[0],hhh1,marker='o',linestyle=' ')
        ax.plot(dff[1],hhh2,marker='s',linestyle=' ')
        ax.plot(dff[2],hhh3,marker='v',linestyle=' ')
        ax.plot(dff[3],hhh4,marker='^',linestyle=' ')
        ax.plot(dff[4],hhh5,marker='+',linestyle=' ')
        ax.plot(dff[5],hhh6,marker='x',linestyle=' ')
        ax.plot(dff[0],yfit11,linestyle='--', label='adsR (m='+str(mx1[0])+')')
        ax.plot(dff[1],yfit22,linestyle='--', label='deh1 (m='+str(mx1[1])+')')
        ax.plot(dff[2],yfit33,linestyle='--', label='deh2 (m='+str(mx1[2])+')')
        ax.plot(dff[3],yfit44,linestyle='--', label='trns (m='+str(mx1[3])+')')
        ax.plot(dff[4],yfit55,linestyle='--', label='desP (m='+str(mx1[4])+')')
        ax.plot(dff[5],yfit66,linestyle='--', label='desH (m='+str(mx1[5])+')')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{H2})$')
        plt.title('RDS for '+title+' Kinetics')
        ax.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(H2_G)_scatterbestfitline.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(dff[0],ppp1,marker='o',linestyle=' ')
        ay.plot(dff[1],ppp2,marker='s',linestyle=' ')
        ay.plot(dff[2],ppp3,marker='v',linestyle=' ')
        ay.plot(dff[3],ppp4,marker='^',linestyle=' ')
        ay.plot(dff[4],ppp5,marker='+',linestyle=' ')
        ay.plot(dff[5],ppp6,marker='x',linestyle=' ')
        ay.plot(dff[0],yfit1,linestyle='--', label='adsR (m='+str(mx[0])+')')
        ay.plot(dff[1],yfit2,linestyle='--', label='deh1 (m='+str(mx[1])+')')
        ay.plot(dff[2],yfit3,linestyle='--', label='deh2 (m='+str(mx[2])+')')
        ay.plot(dff[3],yfit4,linestyle='--', label='trns (m='+str(mx[3])+')')
        ay.plot(dff[4],yfit5,linestyle='--', label='desP (m='+str(mx[4])+')')
        ay.plot(dff[5],yfit6,linestyle='--', label='desH (m='+str(mx[5])+')')
        plt.xlabel('$ln(k_f)$')
        plt.ylabel('$ln(r_{P})$')
        plt.title('RDS for '+title+' Kinetics')
        ay.legend()
        plt.show()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-RDS-k/'+title+'_plot_data_for_RDSk(P_G)_scatterbestfitline.png')

print('COMPLETED')
    
print('CALCULATION COMPLETED')



"""
