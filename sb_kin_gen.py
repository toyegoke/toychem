from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdata
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5
from sb_kin_temp0 import kineticparameter, dydts,solve_odes,solve_odes2,solve_odes3,solve_odes4,solve_odes45,solve_odes46cc
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
# MICRO-KINECTIC & REACTOR MODELS (SEMIBATCH)
#====================================================================================================

def RDS(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowrate, initial_coverage):      # <=============  

    #for specific time range GOOD FOR kinetic simulation  
    
    #for specific time range GOOD FOR RDS
    
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

    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)

    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)

    fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T)     

    gat = [0,0,0,0,0,0,0, 0,0,0];         bat = [1,1,1,1,1,1];         fat = bat;    rr0 = n_R0   

    def RDS_effect2(moleculer_property_function, title_of_molecule, Order_specieID_list, fat, bat, gat): # <=============    
        CCL = []        #RR = []
        RP = []
        RH = []
        T = temp
             
        #xx1,y1= solve_odes46(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)
        xx1,y1= solve_odes46cc(SA,n_R0,v_R0,C_R0,outflow_rate,temp,TIME,dydts, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat, initial_coverage)        

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
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_1amount_plot_.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.plot(xx1, y1[:,6]/y1[0,6], label='$n_R$')
        ay.plot(xx1, y1[:,7]/y1[0,6], label='$n_{H2,in}$')
        #ay.plot(xx1, (y1[:,8]-y1[:,7])/y1[0,6],  label='$n_{H2,out}$')
        ay.plot(xx1, y1[:,8]/y1[0,6], label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ay.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_amount_plot_.png')

        """
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, y1[:,6]/y1[0,6], label='$n_R$')
        ax.semilogx(xx1, y1[:,7]/y1[0,6], label='$n_{H2,in}$')
        ax.semilogx(xx1, (y1[:,8]-y1[:,7])/y1[0,6],  label='$n_{H2,out}$')
        ax.semilogx(xx1, y1[:,8]/y1[0,6], label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_1amount_plot_.png')


        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx1, y1[:,6]/y1[0,6],marker='^',  label='$n_R$')
        ax.semilogx(xx1, y1[:,7]/y1[0,6],marker='x', label='$n_{H2,in}$')
        ax.semilogx(xx1, (y1[:,8]-y1[:,7])/y1[0,6],marker='+',  label='$n_{H2,out}$')
        ax.semilogx(xx1, y1[:,8]/y1[0,6],marker='s', label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_amount_plot_.png')

        fig = plt.figure()
        ay = plt.subplot(111)
        ay.semilogx(xx1, y1[:,0]/S_x,marker='x', label='$\Theta_X$')
        ay.semilogx(xx1, y1[:,1]/S_x,marker='x', label='$\Theta_{RX}$')
        ay.semilogx(xx1, y1[:,2]/S_x,marker='x', label='$\Theta_{UX}$')
        ay.semilogx(xx1, y1[:,3]/S_x,marker='x', label='$\Theta_{HX}$')
        ay.semilogx(xx1, y1[:,4]/S_x,marker='x', label='$\Theta_{VX}$')
        ay.semilogx(xx1, y1[:,5]/S_x,marker='x', label='$\Theta_{PX}$')
        #plt.title(title+' Coverage Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
        ay.legend()  
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_surface_plot_.png')
        """

        fig = plt.figure()
        az = plt.subplot(111)
        az.semilogx(xx1, y1[:,0]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_X$')
        az.semilogx(xx1, y1[:,1]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_{RX}$')
        az.semilogx(xx1, y1[:,2]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_{UX}$')
        az.semilogx(xx1, y1[:,3]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_{HX}$')
        az.semilogx(xx1, y1[:,4]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_{VX}$')
        az.semilogx(xx1, y1[:,5]/(y1[:,0]+y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5]), label='$\Theta_{PX}$')
        #plt.title(title+' Coverage Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Coverage, $ {^{\Theta_i}/_{\Theta_{X,0}}} $')   
        az.legend()  
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_1surface_plot_.png')

        # coverage and concentration definitions 1 (from back) in umol
        n_R = y1[-1,6];  n_H2 = y1[-1,7]; n_P = y1[-1,8]
        Ttim = xx1[-1]; n_Ro = y1[0,6]

        # coverage and concentration definitions 2 (from back) in umol
        n_Rb = y1[-2,6];  n_H2b = y1[-2,7]; n_Pb = y1[-2,8]
        Ttimb = xx1[-2]

        # coverage and concentration definitions 3 (from back) in umol
        utim = int(len(xx1)*0.84); n_Rc = y1[utim,6];  n_H2c = y1[utim,7]; n_Pc = y1[utim,8]
        Ttimc =xx1[utim]

        # coverage definitions 1 (from back) in umol/mol
        x_R = y1[-1,1];  x_U = y1[-1,2];  x_H = y1[-1,3];  x_V = y1[-1,4]; x_P = y1[-1,5]
        Ttim = xx1[-1]; x_Xo = y1[0,0]; x_X = y1[-1,0]

        # coverage definitions 2 (from back) in umol/mol
        x_Rb = y1[-2,1];  x_Ub = y1[-2,2];  x_Hb = y1[-2,3];  x_Vb = y1[-2,4]; x_Pb = y1[-2,5]
        Ttimb = xx1[-2]

        xRR=((x_R - x_Rb)/(Ttim - Ttimb)) # in umol/m2.s 
        xPP=((x_P - x_Pb)/(Ttim - Ttimb)) # in umol/m2.s  
        xHH=((x_H - x_Hb)/(Ttim - Ttimb)) # in umol/m2.s 
        xUU=((x_U - x_Ub)/(Ttim - Ttimb)) # in umol/m2.s 
        xVV=((x_V - x_Vb)/(Ttim - Ttimb)) # in umol/m2.s 

        RRR=((n_R/n_Ro - n_Rb/n_Ro)/(Ttim - Ttimb)) # in 1/s 
        RPP=((n_P/n_Ro - n_Pb/n_Ro)/(Ttim - Ttimb)) 
        RHH=((n_H2/n_Ro - n_H2b/n_Ro)/(Ttim - Ttimb))

        RRR1=((n_R/n_Ro - n_Rc/n_Ro)/(Ttim - Ttimc)) # in 1/s (at 50% time i.e. from -12 to 4 = 12+1.5/16)
        RPP1=((n_P/n_Ro - n_Pc/n_Ro)/(Ttim - Ttimc)) 
        RHH1=((n_H2/n_Ro - n_H2c/n_Ro)/(Ttim - Ttimc))

        RRR2=((n_R/n_Ro - n_Ro/n_Ro)/(Ttim - 0)) # in 1/s (at 50% time i.e. from -12 to 4 = 12+1.5/16)
        RPP2=((n_P/n_Ro - 0/n_Ro)/(Ttim - 0)) 
        RHH2=((n_H2/n_Ro - 0/n_Ro)/(Ttim - 0))

        mk_rate = RPP*n_Ro*(1e-6) # in mol/s
        TOF_mk = mk_rate/(S_x*SA)  # in 1/s
        ES_eff_mk = -R*T*np.log(h*TOF_mk/kB/T)/1000 # in kJ/mol

        mk_rate2 = RPP2*n_Ro*(1e-6) # in mol/s
        TOF_mk2 = mk_rate2/(S_x*SA)  # in 1/s
        ES_eff_mk2 = -R*T*np.log(h*TOF_mk2/kB/T)/1000 # in kJ/mol

        mk_rate1 = outflow_rate/v_R0*(n_P-n_H2)*(1e-6) # in mol/s
        TOF_mk1 = mk_rate1/(S_x*SA)  # in 1/s
        ES_eff_mk1 = -R*T*np.log(h*TOF_mk1/kB/T)/1000 # in kJ/mol

        glen=len(xx1)
        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_amount_data.cvs"
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerow(['Ri_mol/s'+str(T),'Pi_mol/s','Hi_mol/s','Time_s (using average rate=dr/dt)'])
            wr.writerow([ RRR*n_Ro*(1e-6), RPP*n_Ro*(1e-6), RHH*n_Ro*(1e-6) ])
            wr.writerow(['  '])
            #wr.writerow(['Ri_mol/s'+str(T),'Pi_mol/s','Hi_mol/s','Time_s (using average rate=dr/dt from t)'])
            #wr.writerow([ RRR1*n_Ro*(1e-6), RPP1*n_Ro*(1e-6), RHH1*n_Ro*(1e-6) ])
            #wr.writerow(['  '])
            wr.writerow(['Ri_mol/s'+str(T),'Pi_mol/s','Hi_mol/s','Time_s (using average rate=dr/dt from 0)'])
            wr.writerow([ RRR2*n_Ro*(1e-6), RPP2*n_Ro*(1e-6), RHH2*n_Ro*(1e-6) ])
            wr.writerow(['  '])
            wr.writerow(['xRi_mol/m2.s'+str(T),'xUi_mol/m2.s','xHi_mol/m2.s','xVi_mol/m2.s','xPi_mol/m2.s'])
            wr.writerow([ xRR,                xUU,             xHH,             xVV,           xPP ])
            wr.writerow(['  '])
            wr.writerow(['n_H2(mol)=',n_H2*(1e-6), 'x_Xo(mol/m2)=',x_Xo*(1e-6), 'x_X(mol/m2)=',x_X*(1e-6)])
            wr.writerow(['n_Ro(mol)=', n_Ro*(1e-6), '  n_R(mol)=', n_R*(1e-6), '  n_P(mol)=',n_P*(1e-6),'   CONV (X) =', n_P/n_Ro, ' v_R0(m3)=', v_R0])
            wr.writerow(['x_R(mol/m2)=', x_R*(1e-6), 'x_H(mol/m2)=', x_H*(1e-6),'x_U(mol/m2)=', x_U*(1e-6),'x_V(mol/m2)=', x_V*(1e-6),'x_P(mol/m2)=', x_P*(1e-6)])
            wr.writerow(['  '])
            wr.writerow(['n_Ro(mol)=',n_Ro*(1e-6), '  n_R(mol)=',n_R*(1e-6), '  n_P(mol)=',n_P*(1e-6), '  n_H2(mol)=',n_H2*(1e-6),'   CONV (X) =', n_P/n_Ro, ' v_R0(m3)=', v_R0])
            wr.writerow(['  '])
            wr.writerow(['surf conc(mol/m2)','specif surface(m2/g)','initial_mol_R(mol)','mass(g)','flowrate(m3/s)', 'planck const'])
            wr.writerow([S_x, dc, n_R0, mc, outflow_rate, h])        
            wr.writerow(['  '])
            wr.writerow(['bolzmans const','gas const','reactor_vol','exit amount(mol)','inti_conc','cat_area(m2)','amount of site(mol)'])
            wr.writerow([ kB, R, v_R0, n_P*(1e-6), C_R0, SA, S_x*SA])        
            wr.writerow(['  '])
            wr.writerow(['MK rate (mol/s) = dr/dt','TOF (1/s) = mk_rate/(S_x*SA)','ES_eff (kJ/mol) = -RTln(h*TOF/kB/T)'])
            wr.writerow([ mk_rate, TOF_mk, ES_eff_mk ])               
            wr.writerow(['  '])
            wr.writerow(['MK rate (mol/s) = dr/dt from 0','TOF (1/s) = mk_rate/(S_x*SA)','ES_eff (kJ/mol) = -RTln(h*TOF/kB/T)'])
            wr.writerow([ mk_rate2, TOF_mk2, ES_eff_mk2 ])               
            wr.writerow(['  '])
            wr.writerow(['MK rate (mol/s) = vo/vR*nH2','TOF (1/s) = mk_rate/(S_x*SA)','ES_eff (kJ/mol) = -RTln(h*TOF/kB/T)'])
            wr.writerow([ mk_rate1, TOF_mk1, ES_eff_mk1 ])               
            wr.writerow(['  '])
            #wr.writerow(['Ri/Ro'+str(T),'Hi/Ro ','Pi/Ro','Time'])
            wr.writerow(['XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX'])
            wr.writerow(['XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX','XXXXXXXXXXXXXXXXXX'])
            #i=0
            #for i in range(glen):            
            #    wr.writerow([y1[i,6]/y1[0,6], y1[i,7]/y1[0,6], y1[i,8]/y1[0,6], xx1[i]])  

        """
        xx2,y2= solve_odes47(SA,n_R0,v_R0,C_R0,outflow_rate,temp,10**(-4),dydt, G_m, int_id, Order_specieID_list,S_x, fat, bat, gat)
        rr1=y2[:,6]#*1e-6
        hh1=y2[:,7]#*1e-6
        pp1=y2[:,8]#*1e-6

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx2, y2[:,8]/y2[0,6], label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_21amount_plot_.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx2, y2[:,8]/y2[0,6],marker='x', label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_22amount_plot_.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx2, y2[:,6]/y2[0,6],marker='x', label='$n_R$')
        ax.semilogx(xx2, y2[:,8]/y2[0,6],marker='x', label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_23amount_plot_.png')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.semilogx(xx2, y2[:,6]/y2[0,6], label='$n_R$')
        ax.semilogx(xx2, y2[:,8]/y2[0,6], label='$n_P$')
        #plt.title(title+' Production Profile')
        plt.xlabel('Time, $ t/sec $')
        plt.ylabel('Relative Amount, $ {^{n_i}/_{n_{R,0}}} $')
        ax.legend()
        fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/'+title+'_24amount_plot_.png')

        #fig = plt.figure()
        #ax = plt.subplot(111)
        #ax.semilogx(xx2, y2[:,8]/y2[0,6], label='$n_P$_')
        #plt.title('Yield Profile')
        #plt.xlabel('Time, $ t/sec $')
        #plt.ylabel('Yield, $ {^{n_i}/_{n_{R,0}}} $')
        #ax.legend()
        #fig.savefig('/home/toyegoke/ENS_kinetics/Final-gen/COMBINED_25amount_plot_.png')

        print(title+' results')
        print('==================================')
        print('n(P)=',pp1[-1],'mol')
        print('n(H2)=',hh1[-1],'mol')
        print('n(R)=',rr1[-1],'mol')
        print('Y(P)=',pp1[-1]/rr0)
        print('Y(H2)=',hh1[-1]/rr0)
        print('X(R)=',(rr0-rr1[-1])/rr0)
        print('==================================')
      
        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_surface_data.txt"
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerow(['X_'+str(T),'RX  ','UX  ','HX  ','VX  ','PX  ','Time'])
            wr.writerow([y1[1,0],y1[1,1],y1[1,2],y1[1,3],y1[1,4],y1[1,5], xx1])
            wr.writerow([y1[2,0],y1[2,1],y1[2,2],y1[2,3],y1[2,4],y1[2,5], xx1])
            wr.writerow([y1[3,0],y1[3,1],y1[3,2],y1[3,3],y1[3,4],y1[3,5], xx1])
            wr.writerow([y1[4,0],y1[4,1],y1[4,2],y1[4,3],y1[4,4],y1[4,5], xx1])
            wr.writerow([y1[5,0],y1[5,1],y1[5,2],y1[5,3],y1[5,4],y1[5,5], xx1])
            wr.writerow([y1[6,0],y1[6,1],y1[6,2],y1[6,3],y1[6,4],y1[6,5], xx1])
            wr.writerow([y1[7,0],y1[7,1],y1[7,2],y1[7,3],y1[7,4],y1[7,5], xx1])
            wr.writerow([y1[8,0],y1[8,1],y1[8,2],y1[8,3],y1[8,4],y1[8,5], xx1])
            wr.writerow([y1[9,0],y1[9,1],y1[9,2],y1[9,3],y1[9,4],y1[9,5], xx1])
            wr.writerow(['X_'+str(T),'RX  ','UX  ','HX  ','VX  ','PX  ','Time'])
            wr.writerow([y1[-9,0],y1[-9,1],y1[-9,2],y1[-9,3],y1[-9,4],y1[-9,5], xx1])
            wr.writerow([y1[-8,0],y1[-8,1],y1[-8,2],y1[-8,3],y1[-8,4],y1[-8,5], xx1])
            wr.writerow([y1[-7,0],y1[-7,1],y1[-7,2],y1[-7,3],y1[-7,4],y1[-7,5], xx1])
            wr.writerow([y1[-6,0],y1[-6,1],y1[-6,2],y1[-6,3],y1[-6,4],y1[-6,5], xx1])
            wr.writerow([y1[-5,0],y1[-5,1],y1[-5,2],y1[-5,3],y1[-5,4],y1[-5,5], xx1])
            wr.writerow([y1[-4,0],y1[-4,1],y1[-4,2],y1[-4,3],y1[-4,4],y1[-4,5], xx1])
            wr.writerow([y1[-3,0],y1[-3,1],y1[-3,2],y1[-3,3],y1[-3,4],y1[-3,5], xx1])
            wr.writerow([y1[-2,0],y1[-2,1],y1[-2,2],y1[-2,3],y1[-2,4],y1[-2,5], xx1])
            wr.writerow([y1[-1,0],y1[-1,1],y1[-1,2],y1[-1,3],y1[-1,4],y1[-1,5], xx1])

        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_amount_data.txt"
        wr = open(csvfile, "w")
        wr.write('Ri/Ro_'+str(T)+'  Hi/Ro  '+'  Pi/Ro  '+'  Time')
        wr.write(str(y1[i,6]/y1[0,6])+'_'+str(y1[i,7]/y1[0,6])+'_'+str(y1[i,8]/y1[0,6])+'_'+str(xx1[i]))  
        wr.write(str(y1[i,6]/y1[0,6])+'_'+str(y1[i,7]/y1[0,6])+'_'+str(y1[i,8]/y1[0,6])+'_'+str(xx1[i]))
        wr.close()  

        csvfile1 = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_surface_data.txt"
        wr1 = open(csvfile1, "w")
        wr1.write(str('X_'+str(T)+'RX  '+'UX  '+'HX  '+'VX  '+'PX  '+'Time')
        wr1.write(str(y1[-2,0])+'_'+str(y1[-2,1])+'_'+str(y1[-2,2])+'_'+str(y1[-2,3])+'_'+str(y1[-2,4])+'_'+str(y1[-2,5])+'_'+str(xx1))
        wr1.write(str(y1[-1,0])+'_'+str(y1[-1,1])+'_'+str(y1[-1,2])+'_'+str(y1[-1,3])+'_'+str(y1[-1,4])+'_'+str(y1[-1,5])+'_'+str(xx1))
        wr1.close()
 
        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_surface_data.txt"
        wr = open(csvfile, "w")
        wr.write('X_'+str(T),'RX  ','UX  ','HX  ','VX  ','PX  ','Time','\n')
        wr.write(y1[-1,0],y1[-1,1],y1[-1,2],y1[-1,3],y1[-1,4],y1[-1,5], xx1,'\/n')

        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_surface_data.txt"
        wr1 = open(csvfile, "w")
        wr1.write(str('X_'+str(T)+'RX  '+'UX  '+'HX  '+'VX  '+'PX  '+'Time\n')
        wr1.write(str(y1[-2,0])+"_"+str(y1[-2,1])+"_"+str(y1[-2,2])+"_"+str(y1[-2,3])+"_"+str(y1[-2,4])+"_"+str(y1[-2,5])+"_"+str(xx1)+"\n")
        wr1.write(str(y1[-1,0])+"_"+str(y1[-1,1])+"_"+str(y1[-1,2])+"_"+str(y1[-1,3])+"_"+str(y1[-1,4])+"_"+str(y1[-1,5])+"_"+str(xx1)+"\n")
        wr1.close()

        csvfile1 = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_surface_data.txt"
        wr = open(csvfile1, "w")
        wr.write(str('X_'+str(T)+'RX  '+'UX  '+'HX  '+'VX  '+'PX  '+'Time\n')
        wr.write(str(y1[-2,0])+"_"+str(y1[-2,1])+"_"+str(y1[-2,2])+"_"+str(y1[-2,3])+"_"+str(y1[-2,4])+"_"+str(y1[-2,5])+"_"+str(xx1[-2])\n)
        wr.write(str(y1[-1,0])+"_"+str(y1[-1,1])+"_"+str(y1[-1,2])+"_"+str(y1[-1,3])+"_"+str(y1[-1,4])+"_"+str(y1[-1,5])+"_"+str(xx1[-1])\n)
        wr.close()

        wr.write(str(y1[-5,6]/y1[0,6])+"_"+str(y1[-5,7]/y1[0,6])+"_"+str(y1[-5,8]/y1[0,6])+"_"+str(xx1[-5])+"\n")  
        wr.write(str(y1[-4,6]/y1[0,6])+"_"+str(y1[-4,7]/y1[0,6])+"_"+str(y1[-4,8]/y1[0,6])+"_"+str(xx1[-4])+"\n")  
        wr.write(str(y1[-3,6]/y1[0,6])+"_"+str(y1[-3,7]/y1[0,6])+"_"+str(y1[-3,8]/y1[0,6])+"_"+str(xx1[-3])+"\n")  
        wr.write(str(y1[-2,6]/y1[0,6])+"_"+str(y1[-2,7]/y1[0,6])+"_"+str(y1[-2,8]/y1[0,6])+"_"+str(xx1[-2])+"\n")  
        wr.write(str(y1[-1,6]/y1[0,6])+"_"+str(y1[-1,7]/y1[0,6])+"_"+str(y1[-1,8]/y1[0,6])+"_"+str(xx1[-1])+"\n")

        wr.write('        X/X_'+str(T)+'         RX/X         '+'         UX/X        '+'       HX/X        '+'       VX/X       '+'       PX/X       '+'       Time\n')
        wr.write(str(y1[-5,0]/y1[0,0])+"_"+str(y1[-5,1]/y1[0,0])+"_"+str(y1[-5,2]/y1[0,0])+"_"+str(y1[-5,3]/y1[0,0])+"_"+str(y1[-5,4]/y1[0,0])+"_"+str(y1[-5,5]/y1[0,0])+"_"+str(xx1[-5])+"\n")  
        wr.write(str(y1[-4,0]/y1[0,0])+"_"+str(y1[-4,1]/y1[0,0])+"_"+str(y1[-4,2]/y1[0,0])+"_"+str(y1[-4,3]/y1[0,0])+"_"+str(y1[-4,4]/y1[0,0])+"_"+str(y1[-4,5]/y1[0,0])+"_"+str(xx1[-4])+"\n")  
        wr.write(str(y1[-3,0]/y1[0,0])+"_"+str(y1[-3,1]/y1[0,0])+"_"+str(y1[-3,2]/y1[0,0])+"_"+str(y1[-3,3]/y1[0,0])+"_"+str(y1[-3,4]/y1[0,0])+"_"+str(y1[-3,5]/y1[0,0])+"_"+str(xx1[-3])+"\n")  
        wr.write(str(y1[-2,0]/y1[0,0])+"_"+str(y1[-2,1]/y1[0,0])+"_"+str(y1[-2,2]/y1[0,0])+"_"+str(y1[-2,3]/y1[0,0])+"_"+str(y1[-2,4]/y1[0,0])+"_"+str(y1[-2,5]/y1[0,0])+"_"+str(xx1[-2])+"\n")  
        wr.write(str(y1[-1,0]/y1[0,0])+"_"+str(y1[-1,1]/y1[0,0])+"_"+str(y1[-1,2]/y1[0,0])+"_"+str(y1[-1,3]/y1[0,0])+"_"+str(y1[-1,4]/y1[0,0])+"_"+str(y1[-1,5]/y1[0,0])+"_"+str(xx1[-1])+"\n")

        csvfile = "/home/toyegoke/ENS_kinetics/Final-gen/"+title+"_amount_data.txt"
        wr = open(csvfile, "w")
        wr.write('  Ri/Ro_'+str(T)+'    Hi/Ro   '+'    Pi/Ro   '+'   Time\n')
        glen=len(xx1)
        i=0
        for i in range(glen):  
            wr.write(str(y1[i,6]/y1[0,6])+"_"+str(y1[i,7]/y1[0,6])+"_"+str(y1[i,8]/y1[0,6])+"_"+str(xx1[i])+"\n")  

        wr.write('        X/X_'+str(T)+'         RX/X         '+'         UX/X        '+'       HX/X        '+'       VX/X       '+'       PX/X       '+'       Time\n')
        for i in range(glen):         
            wr.write(str(y1[i,0]/y1[0,0])+"_"+str(y1[i,1]/y1[0,0])+"_"+str(y1[i,2]/y1[0,0])+"_"+str(y1[i,3]/y1[0,0])+"_"+str(y1[i,4]/y1[0,0])+"_"+str(y1[i,5]/y1[0,0])+"_"+str(xx1[i])+"\n")  

        wr.close()  
        """ 

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

oldone=[0.001, 0.001, 0, 0.013, 0.985, 0]
#oldone=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
initial_coverage = [oldone, oldone, oldone, oldone,oldone] # [X0, RX1, UX2, HX3, VX4, PX5]  
#old one = [0.001, 0.001, 0, 0.013, 0.985, 0] 

TIM = 3600*10 #* 15 # in sec ==================================================================================3600 sec

pt = [ppt1,ppt2,ppt3,ppt4,ppt5]
titl = [title1,title2,title3,title4,title5]
specie_ID = [engg1,engg2,engg3,engg4,engg5]
#rtime = [10**(TIM),10**(TIM),10**(TIM),10**(TIM),10**(TIM)]
rtime = [TIM,TIM,TIM,TIM,TIM] # in stead old time = 10**(5)] #for 5%/20%   orignal

#rtemp = [530,530,530,530,530]
#rtemp = [540,540,540,540,540]
rtemp = [550,550,550,550,550]
#rtemp = [560,560,560,560,560]
#rtemp = [570,570,570,570,570]

cmass = 750 # 9.9e6 # 1e7 # 1e-3 # in mg ======================================================================== 1e-3 # mg

flowrate = 1 # batch will 0 while semibatch will be 1
ranP = len(titl)
i=0

#"""
for i in range(ranP):
        ppt=pt[i]
        title=titl[i]
        engg=specie_ID[i]
        RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowrate,initial_coverage[i])
        #=====================================================================================================
        #=====================================================================================================
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"""

i=5-1
ppt=pt[i]
title=titl[i]
engg=specie_ID[i]
RDS(ppt,title,engg,rtime[i],rtemp[i],cmass,flowrate,initial_coverage[i])

"""

