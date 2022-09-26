from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdata
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from CR_kin_temp import kineticparameter, dydtc,solve_odes46cc
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

def EaApp(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout,initial_coverage):      # <=============  

    #for specific time range GOOD FOR RDS (for RDS when rate is measured in terms of Rate_x)

    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime # in sec 
    temp = rtemp # in K 
    mc0 = cmass # in mg, milligrams
    C_R0 = 100 # 997.774 # 1 # 1 # mol/m3 =====================================================================================997.74 (original)
    outflow_rate = flowrate * 0.5e-3 #5.5e-12# 1 # 1.5 # in m3/sec ======================================================================14.27*1e-7 (original)
    v_R0 = 0.5 # 14.5e-10 # 1 # 1 # in m3 (45 mL) ==================================================================================45*1e-6 (original)
    n_R0 = v_R0*C_R0 # in mol  

    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K     
    ppt = moleculer_property_function # ppt1, ppt2, where the thermo ppt is computed
    
    #Catalyst properties
    dc = 10 # 0.10 # 350 # 400 #29 # specific surface area of catalyst in m^2/g==================================================================================29 (original)
    mc = mc0/1000 # mass of catalyst in g  ==================================================================mc0/1000
    SA = dc * mc # surface area of catalyst in m^2    
    S_x = 6e-7 #1e-5 #1.8e-7 # 1.9e-9 # 10e-9 # mol/m^2 surface conc ==========================================================18*1e-6 (original)    

#=====================
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    

    def EA_effect(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat):      # <=============  

        tt=[]; xxcat=[]; rrx=[]; uux=[]; hhx=[]; vvx=[]; ppx=[]
        rr1=[]; pp1=[]; hh1=[]; rr2=[]; pp2=[]; hh2=[]
        
        tempt = []         #RR = []
        RP = []; RPP=[]; RP1=[]; RH1=[]; rh=[]; rp=[]; ttim=[]
        RH = []; RHH=[]

        TT = [temp-6,temp-4,temp-2,temp,temp+2,temp+4,temp+6] # USING +10/+8/-20K (exclude last 3 pt)

        tlens = len(TT)
        i = 0
        
        for i in range(tlens):
            
            T=TT[i]
            
            filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-Ea/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
            filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-Ea/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)
        
            G_m, int_id, S_m, H_m, fig, fig2 = ppt(T, filename_link_fomat, filename_link_fomat2)

            fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
            fig2.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'_'+str(int(T))+'K_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)        
            
            tt0, y= solve_odes46cc(SA,n_R0,v_R0,C_R0,outflow_rate,T,TIME,dydtc, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat, initial_coverage)

            #========================================
            xx1=tt0  # time in sec  
            
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.semilogx(xx1, y[:,6]/y[0,6], label='$n_R$')
            ax.semilogx(xx1, y[:,7]/y[0,6], label='$n_{H2,in}$')
            ax.semilogx(xx1, (y[:,8]-y[:,7])/y[0,6],  label='$n_{H2,out}$')
            ax.semilogx(xx1, y[:,8]/y[0,6], label='$n_P$')
            #plt.title(title+' Production Profile')
            plt.xlabel('Time, $ t/sec $')
            plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
            ax.legend()
            fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'_'+str(int(T))+'K_'+'_1amount_plot_.png')
    
            fig = plt.figure()
            ay = plt.subplot(111)
            ay.semilogx(xx1, y[:,0]/S_x, label='$\Theta_X$')
            ay.semilogx(xx1, y[:,1]/S_x, label='$\Theta_{RX}$')
            ay.semilogx(xx1, y[:,2]/S_x, label='$\Theta_{UX}$')
            ay.semilogx(xx1, y[:,3]/S_x, label='$\Theta_{HX}$')
            ay.semilogx(xx1, y[:,4]/S_x, label='$\Theta_{VX}$')
            ay.semilogx(xx1, y[:,5]/S_x, label='$\Theta_{PX}$')
            plt.title(title+' Coverage Profile')
            #plt.xlabel('Time, $ t/sec $')
            plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
            ay.legend()  
            fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'_'+str(int(T))+'K_'+'_1surface_plot_.png')
            #========================================
            
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
        
            rP = SA*1e-6*kf4 * PX / (n_R0/(1e-6)) # in a sec
            rH = SA*1e-6*kf5 * HX**2 / (n_R0/(1e-6))  # in a sec

            #RR.append(rr)
            RP.append(rP)
            RH.append(rH)
            tempt.append(T)

            rP1 = SA*1e-6*(kf4 * PX - kb4*n_P * X)/(n_R0/(1e-6))  # in sec
            rH1 = SA*1e-6*(kf5 * HX**2 - kb5*n_H2 * X**2)/(n_R0/(1e-6))  # in sec

            RP1.append(rP1)
            RH1.append(rH1)   
            
            #RPP.append(-outflow_rate/v_R0*(0-n_P))  #rateP= -outflow_rate/v_R0*(0-n_P)
            #RHH.append(-outflow_rate/v_R0*(0-n_H2))  # rateH2= -outflow_rate/v_R0*(0-n_H2)
            RPP.append(-outflow_rate/v_R0*(0-n_P*(1e-6))/n_R0)  #rateP= -outflow_rate/v_R0*(0-n_P)
            RHH.append(-outflow_rate/v_R0*(0-n_H2*(1e-6))/n_R0)  # rateH2= -outflow_rate/v_R0*(0-n_H2)

            #RPP.append((n_P/n_Ro - n_Pb/n_Ro)/(Ttim - Ttimb)) 
            #RHH.append((n_H2/n_Ro - n_H2b/n_Ro)/(Ttim - Ttimb))
            
            # case 1 when gas specie is in n/V, mol/m3.
            rf1 = kf1 * n_R/v_R0 * X
            rb1 = kb1 * RX  
            r1 = SA*1e-6*(rf1)/n_R0#-(rb1)/n_R0
            
            rf2 = kf2 * RX * X  
            rb2 = kb2 * UX * HX  
            r2 = SA*1e-6*(rf2)/n_R0#-(rb2)/n_R0
            
            rf3 = kf3 * UX * X   
            rb3 = kb3 * PX * HX  
            r3 = SA*1e-6*(rf3)/n_R0#-(rb3)/n_R0

            rf6 = kf6 * VX  
            rb6 = kb6 * PX
            r4 = SA*1e-6*(rf6)/n_R0#-(rb6)/n_R0
            
            rf4 = kf4 * PX  
            rb4 = kb4 * n_P/v_R0 * X
            r5 = SA*1e-6*(rf4)/n_R0#-(rb4)/n_R0
            
            rf5 = kf5 * HX**2   
            rb5 = kb5 * n_H2/v_R0 * (X**2)
            r6 = SA*1e-6*(rf5)/n_R0#-(rb5)/n_R0
            
            Rxn=[r1,r2,r3,r4,r5,r6]

            ef1=TSRg-(Rg+Xg)
            eb1=TSRg-(RXg)
            
            ef2=(TS1g+Xg)-(RXg+Xg)
            eb2=(TS1g+Xg)-(UXg+HXg)
            
            ef3=(TS2g+Xg)-(UXg+Xg)
            eb3=(TS2g+Xg)-(VXg+HXg)
            
            ef4=PXg-VXg
            eb4=PXg-PXg
            
            ef5=TSPg-(Pg+Xg)
            eb5=TSPg-PXg
            
            ef6=TSH2g-(H2g+2*Xg)
            eb6=TSH2g-(2*HXg)
            
            Efh=[ef1,ef2,ef3,ef4,ef5,ef6]
            Efp=[ef1,ef2,ef3,ef4,ef5,0]
            Ebh=[eb1,eb2,eb3,eb4,eb5,eb6]          
            Ebp=[eb1,eb2,eb3,eb4,eb5,0]            
            
            Efhmax=max(Efh)
            Efpmax=max(Efp)
            Ebhmax=max(Ebh)
            Ebpmax=max(Ebp)
            
            Ehindex=Efh.index(Efhmax)
            Epindex=Efp.index(Efpmax)

            rhmin = abs(Rxn[Ehindex])  # selectivity for H2
            rpmin = abs(Rxn[Epindex])  # selectivity for P

            rh.append(rhmin)
            rp.append(rpmin)
            
        #print(tempt, rh, rp, loge(rh), loge(rp))
        #print(tempt, RH, RP, loge(RH), loge(RP))
        #print(PX, SA, kf4, n_R0)

        #return tempt, RH, RP, loge(RH), loge(RP)     # desorp rate
        return tempt, RHH, RPP, loge(RHH), loge(RPP) # production profile amount as a rate
        #return tempt, RH1, RP1, loge(RH1), loge(RP1)  # desorp-adsorp rate   
        #return tempt, rh, rp, loge(rh), loge(rp)  # using the step with min kinetic rate   

    def gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list):
       
       title=title_of_molecule
       gat=[0,0,0,0,0,0,0, 0,0,0]
       bat=[1,1,1,1,1,1]
       fat=bat
       ppt = moleculer_property_function # ppt1, ppt2,
       
       tempt, RH, RP, lnRH, lnRP  = EA_effect(ppt, title0, Order_specieID_list, fat, bat, gat)
       slope1, intercept1 = np.polyfit(tempt, lnRH, 1)
       slope2, intercept2 = np.polyfit(tempt, lnRP, 1)
       slope11, intercept11 = np.polyfit(tempt, RH, 1)
       slope22, intercept22 = np.polyfit(tempt, RP, 1)       
       AppEa1 = R*temp**2*slope1/1000
       AppEa2 = R*temp**2*slope2/1000
       AppEa11 = R*temp**2*slope11/1000
       AppEa22 = R*temp**2*slope22/1000

       templen=len(tempt) 
       i=0 
       lnRHfit=[]
       lnRPfit=[]
       RHfit=[]
       RPfit=[]
       for i in range(templen):
           aaa1 = slope1*float(tempt[i]) + intercept1
           aaa2 = slope2*float(tempt[i]) + intercept2
           aaa11 = slope11*float(tempt[i]) + intercept11
           aaa22 = slope22*float(tempt[i]) + intercept22
           lnRHfit.append(aaa1)
           lnRPfit.append(aaa2)
           RHfit.append(aaa11)
           RPfit.append(aaa22)
       """
       # for HYDROGEN
       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, lnRH,marker='o', label='Data Point')
       ax.plot(tempt, lnRHfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa1))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,DP,LF,lnR).png')
       
       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, RH,marker='o', label='Data Point')
       ax.plot(tempt, RH)
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp ')
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,DP,R).png')
       """

       """
       fig = plt.figure()
       ax = plt.subplot(111)
       ax.plot(tempt, lnRHfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa1))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,LF,lnR).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, lnRH,marker='o', label='Data Point')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa1))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,DP,lnR).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, RH,marker='o', label='Data Point')
       ax.plot(tempt, RHfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa11))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,DP,LF,R).png')


       fig = plt.figure()
       ax = plt.subplot(111)
       ax.plot(tempt, RHfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa11))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(H2,LF,R).png')
       """

        # for PRODUCT
       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, lnRP,marker='o', label='Data Point')
       ax.plot(tempt, lnRPfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('At '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa2))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,DP,LF,lnR).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, RP,marker='o', label='Data Point')
       ax.plot(tempt, RP)
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('At '+str(temp)+'K'+' Eapp ')
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,DP,R).png')

       """
       fig = plt.figure()
       ax = plt.subplot(111)
       ax.plot(tempt, lnRPfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa2))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,LF,lnR).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, lnRP,marker='o', label='Data Point')
       plt.xlabel('$T / K$')
       plt.ylabel('$ln(Rate)$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa2))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,DP,lnR).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, RP,marker='o', label='Data Point')
       ax.plot(tempt, RPfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa22))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,DP,LF,R).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.scatter(tempt, RP,marker='o', label='Data Point')
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa22))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,DP,R).png')

       fig = plt.figure()
       ax = plt.subplot(111)
       ax.plot(tempt, RPfit, label='Linear Fit')
       plt.xlabel('$T / K$')
       plt.ylabel('$Rate$')
       plt.title('For '+title+' at '+str(temp)+'K'+' Eapp = %f kJ/mol' % (AppEa22))
       ax.legend()
       #plt.show()
       fig.savefig('/home/toyegoke/ENS_kinetics/Final-Ea/'+title+'__Ea_plot_for_(P,LF,R).png')
       """

       return  AppEa1, AppEa2

    return gat_RDS2(title_of_molecule, moleculer_property_function, Order_specieID_list)


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
#stepRDS = ['3f','3f','4f','3f','5f']

cmass = 750 #1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

rtemp = [550,550,550,550,550] # 550K
#rtemp = [530,530,530,530,530] # 530K
#rtemp = [570,570,570,570,570] # 570K

tim = 3600*10 # since the values is constant for over 10**5
rtime = [tim,tim,tim,tim,tim] #for 5%/20%

oldone=[0.001, 0.001, 0, 0.013, 0.985, 0]
#oldone=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
initial_coverage = [oldone, oldone, oldone, oldone,oldone] # [X0, RX1, UX2, HX3, VX4, PX5]  

"""
ranP = len(titl)
i=0

for i in range(ranP):
    ppt=pt[i]
    title=titl[i]
    engg=specie_ID[i]
    #RDS=stepRDS[i]
    AppEa1, AppEa2 = EaApp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout,initial_coverage[i])
    AppEa1
    AppEa2

print('CALCULATION COMPLETED')

"""

i=1-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
#RDS=stepRDS[i]
AppEa1, AppEa2 = EaApp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout,initial_coverage[i])
AppEa1
AppEa2

i=3-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
#RDS=stepRDS[i]
AppEa1, AppEa2 = EaApp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout,initial_coverage[i])
AppEa1
AppEa2

i=4-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
#RDS=stepRDS[i]
AppEa1, AppEa2 = EaApp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout,initial_coverage[i])
AppEa1
AppEa2

i=5-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
#RDS=stepRDS[i]
AppEa1, AppEa2 = EaApp(ppt,title,engg,rtime[i],rtemp[i],cmass,flowout,initial_coverage[i])
AppEa1
AppEa2

print('CALCULATION COMPLETED')

#"""
     
