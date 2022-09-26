#from gn_ppt import G_m, int_id, S_m, H_m
from gn_thermop import specie_index
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
from numpy import round
import csv

# COMPUTATOR FOR THE ESTIMATION OF TEMP EFFECT ON THE KINETICS OF A REACTION


#=====================================================================================================
#=====================================================================================================
# RATE CONSTANTS CALCULATIONS (BATCH)
#=====================================================================================================
#=====================================================================================================

class kineticparameter:
    def ksu(G_ts_state,G_reacting_species,saturation_surface_conc,gas_phase_specie_mole,surface_specie_mole,T):
        """
        SURFACE REACTION CONSTANT 
        using G, Co(stand. gas conc.), n_gas(specie), n_surf(specie)
        """
        kB=const.Boltzmann # J/K
        N=const.Avogadro # mol^-1
        h=const.Planck # J.s
        R=N*kB #const.gas_constant in J/mol/K
        p=1.01325e5 # Pa
        # standard_gas_conc=p/R/298.15
        # conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        conc_correction = p**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        G_act=max(G_ts_state-G_reacting_species,0) 
        k_rxn=conc_correction*(kB*T/h)*np.exp(-G_act/R/T)
        #print('G_act, k_rxn =',G_act, k_rxn )
        print('k_su(A) =====',conc_correction*(kB*T/h),'(E)=====',G_act)
        return k_rxn

    def kad(Sgas,st_gas,Hgas,dHt_gas,Hads,Hcat,Ggas,Gads,Gcat,saturation_surface_conc,gas_phase_specie_mole,surface_specie_mole,T):
        """
        ADSORPTION RATE CONSTANT calculation
        A is catalyst surface area in m^2
        vg is volume of gas-phase specie in m^3
        Sgas for gas entropy while st_gas is gas trans. for gas entropy
        Hgas for gas enthalpy while dHt_gas is gas trans. for gas enthalpy
        Hcat is enthalpy for free catalyst surface
        Hads is enthalpy for cat. with adsorbed specie
        T is temperature
        """
        K=1 # transmission coefficient factor
        p=1.01325e5 # standard pressure in Pa
        kB=const.Boltzmann
        N=const.Avogadro
        h=const.Planck
        R=N*kB #const.gas_constant
        # standard_gas_conc=p/R/298.15
        # conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        conc_correction = p**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        Sts=Sgas-(1/3*st_gas)
        S_act=Sts-Sgas   # Scat_ts - Scat_free = 0
        Hts=Hgas-(1/3*dHt_gas)
        H_act=Hts-Hgas   # Hcat_ts - Hcat_free = 0
        print('Sts=',Sts,'Hts=',Hts)
        G_act0=H_act-T*S_act
        G_ads=Gads-(Ggas+Gcat)
        G_act=max(0,G_ads,G_act0)
        k_ad=conc_correction*K*(kB*T/h)*np.exp(-G_act/R/T)
        print('S_act=',S_act,'G_act=',G_act,'G_ads=',G_ads,'k_ad =',k_ad)
        print('k_ad (A) = ',conc_correction*K*(kB*T/h),'(E) = ',G_act)
        return k_ad

    def kad2(Sgas,st_gas,Hgas,dHt_gas,Hads,Hcat,Ggas,Gads,Gcat,saturation_surface_conc,gas_phase_specie_mole,surface_specie_mole,T):
        """
        ADSORPTION RATE CONSTANT calculation
        A is catalyst surface area in m^2
        vg is volume of gas-phase specie in m^3
        Sgas for gas entropy while st_gas is gas trans. for gas entropy
        Hgas for gas enthalpy while dHt_gas is gas trans. for gas enthalpy
        Hcat is enthalpy for free catalyst surface
        Hads is enthalpy for cat. with adsorbed specie
        T is temperature
        """
        K=1 # transmission coefficient factor
        p=1.01325e5 # standard pressure in Pa
        kB=const.Boltzmann
        N=const.Avogadro
        h=const.Planck
        R=N*kB #const.gas_constant
        # standard_gas_conc=p/R/298.15
        # conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        conc_correction = p**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        Sts=Sgas-(1/3*st_gas)
        S_act=Sts-Sgas   # Scat_ts - Scat_free = 0
        Hts=Hgas-(1/3*dHt_gas)
        H_act=Hts-Hgas   # Hcat_ts - Hcat_free = 0
        print('Sts=',Sts,'Hts=',Hts)
        H_ads=Hads-(Hgas+Hcat)
        H_act=max(0,H_ads,H_act)
        k_ad=conc_correction*K*(kB*T/h)*np.exp(S_act/R)*np.exp(-H_act/R/T)
        print('S_act=',S_act,'H_act',H_act,'H_ads=',H_ads,'k_ad=',k_ad)
        print('k_ad (A) = ',conc_correction*K*(kB*T/h)*np.exp(S_act/R),'(E) = ',H_act)
        return k_ad

    def kad3(S_ReactSpecies,H_ReactSpecies,G_ReactSpecies,  S_TS,H_TS,G_TS,  S_ads,H_ads,G_ads,  saturation_surface_conc,gas_phase_specie_mole,surface_specie_mole,T):
        """
        ADSORPTION RATE CONSTANT calculation
        A is catalyst surface area in m^2
        vg is volume of gas-phase specie in m^3
        Sgas for gas entropy while st_gas is gas trans. for gas entropy
        Hgas for gas enthalpy while dHt_gas is gas trans. for gas enthalpy
        Hcat is enthalpy for free catalyst surface
        Hads is enthalpy for cat. with adsorbed specie
        T is temperature
        """
        K=1 # transmission coefficient factor
        p=1.01325e5 # standard pressure in Pa
        kB=const.Boltzmann
        N=const.Avogadro
        h=const.Planck
        R=N*kB #const.gas_constant
        # standard_gas_conc=p/R/298.15
        # conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        conc_correction = p**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        dH_ads=H_ads - H_ReactSpecies 
        dH_act=H_TS - H_ReactSpecies 
        H_act=max(0, dH_ads, dH_act)
        dS_ads=S_ads - S_ReactSpecies 
        dS_act=S_TS - S_ReactSpecies 
        S_act=max(0, dS_ads, dS_act)
        k_ad=conc_correction*K*(kB*T/h)*np.exp(S_act/R)*np.exp(-H_act/R/T)
        print('S_act=',S_act,'H_act',H_act,'H_ads=',H_ads,'k_ad=',k_ad)
        print('k_ad (A) = ',conc_correction*K*(kB*T/h)*np.exp(S_act/R),'(E) = ',H_act)
        return k_ad        

    def kad4(S_ReactSpecies,H_ReactSpecies,G_ReactSpecies,  S_TS,H_TS,G_TS,  S_ads,H_ads,G_ads,  saturation_surface_conc,gas_phase_specie_mole,surface_specie_mole,T):
        """
        ADSORPTION RATE CONSTANT calculation
        A is catalyst surface area in m^2
        vg is volume of gas-phase specie in m^3
        Sgas for gas entropy while st_gas is gas trans. for gas entropy
        Hgas for gas enthalpy while dHt_gas is gas trans. for gas enthalpy
        Hcat is enthalpy for free catalyst surface
        Hads is enthalpy for cat. with adsorbed specie
        T is temperature
        """
        K=1 # transmission coefficient factor
        p=1.01325e5 # standard pressure in Pa
        kB=const.Boltzmann
        N=const.Avogadro
        h=const.Planck
        R=N*kB #const.gas_constant
        # standard_gas_conc=p/R/298.15
        # conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        conc_correction = p**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        dG_ads=G_ads - G_ReactSpecies 
        dG_act=G_TS - G_ReactSpecies 
        G_act=max(0, dG_ads, dG_act)

        k_ad=conc_correction*K*(kB*T/h)*np.exp(-G_act/R/T)
        print('G_act',G_act,'G_ads=',G_ads,'k_ad=',k_ad)
        print('k_ad (A) = ',conc_correction*K*(kB*T/h),'(E) = ',G_act)
        return k_ad     

    def keq(Greact,Gprod,T):
        """
        EQUILIBRIUM CONSTANT calculation
        A is catalyst surface area in m^2
        vg is volume of gas-phase specie in m^3
        Ggas for gas free energy 
        Gcat is free energy for free catalyst surface
        Gads is free energy for cat. with adsorbed specie
        T is temperature
        """
        kB=const.Boltzmann
        N=const.Avogadro
        R=N*kB
        G_des=Gprod-Greact
        G_rxn=G_des
        k_equ=np.exp(-G_rxn/R/T) 
        #print('G_rxn, k_eq =',G_rxn, k_equ)
        print('G_rxn = ',G_rxn)
        return k_equ


#=====================================================================================================
#=====================================================================================================
# ODE MODEL SOLVER
#=====================================================================================================
#=====================================================================================================
def solve_odes46b0(SA, n_R0,v_R0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
    """
    Time-integrate chemo-kinetic system
    mc0 (gram)e.g.25, C_R0 (mol/m3)e.g.950, outflow_rate(m3/s)e.g. 5e-7, 
    Temp(K), time(sec), dydt, Energylist, intID_list, Order_specieID_list
    fat=kf_change,bat=kb_change,gat=G_change 
    """
    #fat=[1,1,1,1,1,1]; bat=[1,1,1,1,1,1]; gat=[0,0,0,0,0,0,0]  
    #that is the arguments.....k1,k2,k3,k4,k5,k6..........rx,ux,hx,vx,px,ts1,ts2 

    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML

    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T) i.e. variable
    v_initial = v_R0  # m3, volume the reactor/batch when it is independent on temp (T) i.e. fixed
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate*0.1, v_batch, n_in, SA]

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0*0.001, S_x0*0.001, S_x0*0.00, S_x0*0.013, S_x0*0.985, S_x0*0,   n_R0,   0,    0]  
    #            X,        RX,     UX,           HX,      Vx,          PX,       R,    H,    P
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydtb).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model, fat, bat, gat])

    # integrate on a logaritmic scale
    xx = np.linspace(0, t1, 999) #3990000)
    #xx = np.linspace(-2.0,np.log10(t1), int(t1*18))  #9999900)  #int(np.log10(((t1*2000)+24)*21000))+31000) #900 #800 #600) #int((np.log10(t1-200)+24.0)*100)) #(-14, np.log10(t1), int((np.log10(t1-500)+100)*200)) #int(t1*200))  #t1), 600)#
    
    tt = []
    time = []
    yy = []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        
        if (current[-3]/n_R0)<0.80:   # 0.00001(100%) / 0.80(20%) / 0.95(5%) for time selection
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            break
        
        if tnew>(10**5):
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew,' secs')
            break
        
    print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
    return tt, np.matrix(yy)   

def solve_odes46b1(SA, n_R0,v_R0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
    """
    Time-integrate chemo-kinetic system
    mc0 (gram)e.g.25, C_R0 (mol/m3)e.g.950, outflow_rate(m3/s)e.g. 5e-7, 
    Temp(K), time(sec), dydt, Energylist, intID_list, Order_specieID_list
    fat=kf_change,bat=kb_change,gat=G_change 
    """
    #fat=[1,1,1,1,1,1]; bat=[1,1,1,1,1,1]; gat=[0,0,0,0,0,0,0]  
    #that is the arguments.....k1,k2,k3,k4,k5,k6..........rx,ux,hx,vx,px,ts1,ts2 

    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML

    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T) i.e. variable
    v_initial = v_R0  # m3, volume the reactor/batch when it is independent on temp (T) i.e. fixed
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate*0.1, v_batch, n_in, SA]

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0*1, S_x0*0.000, S_x0*0.000, S_x0*0.000, S_x0*0.000, S_x0*0,   n_R0,   0,    0]  
    #         X,        RX,     UX,           HX,      Vx,          PX,       R,    H,    P
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydtb).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model, fat, bat, gat])

    # integrate on a logaritmic scale
    xx = np.linspace(0, t1, 2)
    #xx = np.linspace(-2.0,np.log10(t1), int(t1*18))  #9999900)  #int(np.log10(((t1*2000)+24)*21000))+31000) #900 #800 #600) #int((np.log10(t1-200)+24.0)*100)) #(-14, np.log10(t1), int((np.log10(t1-500)+100)*200)) #int(t1*200))  #t1), 600)#
    
    tt = []
    time = []
    yy = []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        
        if (current[-3]/n_R0)<0.0050:   # 0.00001(100%) / 0.80(20%) / 0.95(5%) for time selection
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            break
        """
        if tnew>(10**6):
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew,' secs')
            break """
        
    print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
    return tt, np.matrix(yy)    
    
def solve_odes46b(SA, n_R0,v_R0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
    """
    Time-integrate chemo-kinetic system
    mc0 (gram)e.g.25, C_R0 (mol/m3)e.g.950, outflow_rate(m3/s)e.g. 5e-7, 
    Temp(K), time(sec), dydt, Energylist, intID_list, Order_specieID_list
    fat=kf_change,bat=kb_change,gat=G_change 
    """
    #fat=[1,1,1,1,1,1]; bat=[1,1,1,1,1,1]; gat=[0,0,0,0,0,0,0]  
    #that is the arguments.....k1,k2,k3,k4,k5,k6..........rx,ux,hx,vx,px,ts1,ts2 

    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML

    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T) i.e. variable
    v_initial = v_R0  # m3, volume the reactor/batch when it is independent on temp (T) i.e. fixed
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate*0.1, v_batch, n_in, SA]

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0*1, S_x0*0.000, S_x0*0.000, S_x0*0.000, S_x0*0.000, S_x0*0,   n_R0,   0,    0]  
    #         X,        RX,     UX,           HX,      Vx,          PX,       R,    H,    P
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydtb).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model, fat, bat, gat])

    # integrate on a logaritmic scale
    #xx = np.linspace(0, t1, 990000)
    xx = np.linspace(-2.0,np.log10(t1), int(t1*18))  #9999900)  #int(np.log10(((t1*2000)+24)*21000))+31000) #900 #800 #600) #int((np.log10(t1-200)+24.0)*100)) #(-14, np.log10(t1), int((np.log10(t1-500)+100)*200)) #int(t1*200))  #t1), 600)#
    
    tt = []
    time = []
    yy = []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        
        if (current[-3]/n_R0)<0.00500:   # 0.00001(100%) / 0.80(20%) / 0.95(5%) for time selection
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            break
        
        """
        if tnew>(10**6):
            time.append(tnew/60/60)
            print('the reaction time (t_max) was ',tnew,' secs')
            break """
        
    print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
    return tt, np.matrix(yy)   
#=====================================================================================================

#=====================================================================================================
#SEMI-BATCH MODEL
#=====================================================================================================
def dydtb(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
    G_m =  params[1]; energylist = G_m     # energy list
    int_id = params[2];     # list of specie IDs  with string  elements

    rxn_specie_id = params[3]  # list of selected specie IDs in a special order with string  elements
    x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
    vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; 
    ts2=rxn_specie_id[10]; 
    ts3=rxn_specie_id[11]; 
    tsr=rxn_specie_id[12]; 
    tsp=rxn_specie_id[13]; 
    tsh2=rxn_specie_id[14] 

    reactor_model = params[4]; outflow_rate=reactor_model[0]; v_batch=reactor_model[1]; n_in=reactor_model[2]; SA=reactor_model[3]; v_R0=v_batch
    fat = params[5]; fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
    bat = params[6]; bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
    gat = params[7]; gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]; gat8=gat[7]; gat9=gat[8]; gat10=gat[9]; gat11=gat[10]
    dydt = np.zeros(9)

    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    # Co = float(p_o/R/298.15) # in mol/m^3 (standard conc)
    # Co = p_o # in Pa (standard pressure)
    
    # reaction mechanism steps
    """
    2R + 4X <=TSR=> RX + 3X
    RX + 3X <=TS1=> UX + HX + 2X
    UX + HX + 2X <=TS2=> VX + 2HX + X
    VX + 2HX + X <=TS3=> PX + 3HX
    PX + 3HX <=TSP=> P + 3HX + X
    P + 3HX + X <=TSH=> P + 3H2 + 4X
    or
    2R + X <=> RX
    RX + X <=> UX + HX
    UX + X <=> VX + HX
    VX + X <=> PX + HX
    PX <=> P + X
    HX <=> H2 + X
    """

    # coverage and concentration definitions
    X  = y[0]    
    RX = y[1]
    UX = y[2]
    HX = y[3]
    VX = y[4]
    PX = y[5]
    n_R = y[6]
    n_H2 = y[7]
    n_P = y[8]

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

    # Rate of elementary reaction steps mol/m2/sec
    """
    # case 1 when gas specie is in n, mol.  

    rf1 = kf1 * n_R * X
    rb1 = kb1 * RX  
    rf2 = kf2 * RX * X  
    rb2 = kb2 * UX * HX  
    rf3 = kf3 * UX * X   
    rb3 = kb3 * PX * HX  
    rf6 = kf6 * VX  
    rb6 = kb6 * PX
    rf4 = kf4 * PX  
    rb4 = kb4 * (n_P**0.5) * X
    rf5 = kf5 * HX**3   
    rb5 = kb5 * (n_H2**1.5) * (X**3)    
    """
    # case 1 when gas specie is in P = nRT/V, Pa.

    pR = n_R*R*T/v_R0
    pP = n_P*R*T/v_R0
    pH2 = n_H2*R*T/v_R0

    rf1 = kf1 * (pR**2) * X
    rb1 = kb1 * RX  
    rf2 = kf2 * RX * X  
    rb2 = kb2 * UX * HX  
    rf3 = kf3 * UX * X   
    rb3 = kb3 * PX * HX  
    rf4 = kf4 * PX  
    rb4 = kb4 * (pP**1) * X
    rf5 = kf5 * (HX**1)   
    rb5 = kb5 * (pH2**1) * (X**1)
    rf6 = kf6 * VX  
    rb6 = kb6 * PX

    # ODEs for rate of change in surface specie coverages
    dydt[0] = 1*(-rf1 + rb1 - rf2 + rb2 - rf3 + rb3 + rf4 - rb4  + 3*rf5 - 3*rb5) # X in mol/m2/sec
    dydt[1] = 1*(rf1 - rb1 - rf2 + rb2) # RX in mol/m2/sec
    dydt[2] = 1*(rf2 - rb2 - rf3 + rb3) # UX in mol/m2/sec
    dydt[3] = 1*(rf2 - rb2 + rf3 - rb3 - 3*rf5 + 3*rb5) # HX in mol/m2/sec
    dydt[4] = 1*(-rf6 + rb6) # VX in mol/m2/sec
    dydt[5] = 1*(rf3 - rb3 - rf4 + rb4 + rf6 - rb6) # PX in mol/m2/sec

    # ODEs for rate of change in gas-phase concentration BATCH REACTOR
    dydt[6] = SA*1e-6*(-rf1 + rb1) # R in mol/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5)  # H2 in mol/sec
    dydt[8] = SA*1e-6*(rf4 - rb4)  # P in mol/sec

    # display of rate(s) for the steps
    print('rf1=',rf1)
    print('rb1=',rb1) 
    print(rf1-rb1<0,rf1-rb1)   
    print('rf2=',rf2)    
    print('rb2=',rb2)
    print(rf2-rb2>0,rf2-rb2)
    print('rf3=',rf3)
    print('rb3=',rb3)      
    print(rf3-rb3<0,rf3-rb3)
    print('rf6=',rf6)    
    print('rb6=',rb6)
    print(rf6-rb5<6,rf6-rb6)
    print('rf4=',rf4)    
    print('rb4=',rb4) 
    print(rf4-rb4>0,rf4-rb4)
    print('rf5=',rf5)    
    print('rb5=',rb5)
    print(rf5-rb5<0,rf5-rb5)
    print()  

    # display of rate constant(s) for the steps
    print('kf1=',kf1)   
    print('kb1=',kb1)
    print(kf1-kb1<0,kf1-kb1)     
    print('kf2=',kf2)    
    print('kb2=',kb2)
    print(kf2-kb2>0,kf2-kb2)
    print('kf3=',kf3)         
    print('kb3=',kb3)
    print(kf3-kb3<0,kf3-kb3)    
    print('kf6=',kf6)         
    print('kb6=',kb6)
    print(kf6-kb6<0,kf6-kb6)    
    print('kf4=',kf4)    
    print('kb4=',kb4)
    print(kf4-kb4>0,kf4-kb4) 
    print('kf5=',kf5)    
    print('kb5=',kb5)
    print(kf5-kb5<0,kf5-kb5)

    # display of concentration(s) and coverage(s) for the intermediate & gas phase species 
    print('n_R=',n_R,'n_H2=',n_H2, 'n_P=',n_P, 'n_sum=',n_R+n_H2+n_P, 'n_sum11=',n_R+0.5*n_H2+0.5*n_P)
    print('RX=',RX/(X+RX+HX+UX+VX+PX), 'UX=',UX/(X+RX+HX+UX+VX+PX), 'HX=',HX/(X+RX+HX+UX+VX+PX), 'VX=',VX/(X+RX+HX+UX+VX+PX),'PX=',PX/(X+RX+HX+UX+VX+PX),'X=',X/(X+RX+HX+UX+VX+PX),'Xsum=',X+RX+HX+UX+VX+PX)
    print()
    print('DYDT8=',dydt[8])
    print('DYDT8=',dydt[8])
    print('DYDT8=',dydt[8])
    return dydt
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
