#from ppt import G_m, int_id, S_m, H_m
from thermop import specie_index
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
from numpy import round

# COMPUTATOR FOR THE ESTIMATION OF TEMP EFFECT ON THE KINETICS OF A REACTION


#=====================================================================================================
#=====================================================================================================
# RATE CONSTANTS CALCULATIONS
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
        standard_gas_conc=p/R/298.15
        conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
        G_act=G_ts_state-G_reacting_species 
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
        standard_gas_conc=p/R/298.15
        conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
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
        standard_gas_conc=p/R/298.15
        conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
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
        standard_gas_conc=p/R/298.15
        conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
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
        standard_gas_conc=p/R/298.15
        conc_correction=standard_gas_conc**(-gas_phase_specie_mole)*saturation_surface_conc**(1-surface_specie_mole)
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

def solve_odes(T,tim,dydt,Gm,intID,specieID_list, S_x):
    """
    Time-integrate chemo-kinetic system
    Temp, time, dydt, Energylist, intID_list, Order_specieID_list
    """
    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    #Feed parameters
    C_R0 = 950 # in mol/m3 (0.95 M)
    v_R0 = 45*1e-6 # in m3 (45 mL)
    n_R0 = v_R0*C_R0 # in mol    
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = 25/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0
    outflow_rate = 5e-7 # m3/sec, volumetric flowrate for the H2 withdrawal 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim*60*60 # total integration time
    
    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=500, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), int((np.log10(t1) + 0.12)*10))
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
    return tt, yy[0], yy[1], yy[2], yy[3], yy[4], yy[5], yy[6], yy[7], yy[8]     
#=====================================================================================================
def solve_odes2(T,tim,dydt,Gm,intID,specieID_list, S_x):
    """
    Time-integrate chemo-kinetic system
    Temp, time, dydt, Energylist, intID_list, Order_specieID_list
    """
    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    #Feed parameters
    C_R0 = 950 # in mol/m3 (0.95 M)
    v_R0 = 45*1e-6 # in m3 (45 mL)
    n_R0 = v_R0*C_R0 # in mol    
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = 25/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0
    outflow_rate = 5e-7 # m3/sec, volumetric flowrate for the H2 withdrawal 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=500, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), 5)#int((np.log10(t1) + 12)*10e2))
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
    return tt, np.matrix(yy)     
#=====================================================================================================
def solve_odes3(mc0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x):
    """
    Time-integrate chemo-kinetic system
    mc0 (milligram)e.g.25, C_R0 (mol/m3)e.g.950, outflow_rate(x 1e-7 m3/s)e.g. 5, 
    Temp(K), time(sec), dydt, Energylist, intID_list, Order_specieID_list
    """
    p=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam ni
    R=N*kB # const.gas_constant in J/mol/K
    #Feed parameters
    v_R0 = 45*1e-6 # in m3 (45 mL)
    n_R0 = v_R0*C_R0 # in mol    
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=500, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), 2)#int((np.log10(t1) + 12)*10e2))
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
    return tt, np.matrix(yy)     

#=====================================================================================================

#=====================================================================================================
#=====================================================================================================
def solve_odes4(n_R0,v_R0,mc0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
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
  
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model, fat, bat, gat])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), 2)#
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
    return tt, np.matrix(yy)   


#=====================================================================================================
#=====================================================================================================
def solve_odes40(SA, n_R0,v_R0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
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
  
    #Catalyst properties
    #dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    #mc = mc0/1000 # mass of catalyst in g
    #SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, Gm, intID, specieID_list, reactor_model, fat, bat, gat])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), 2)#
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
    return tt, np.matrix(yy)   


#=====================================================================================================

def dydt(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
    G_m =  params[1]; energylist = G_m     # energy list
    int_id = params[2];     # list of specie IDs  with string  elements
    rxn_specie_id = params[3]     # list of selected specie IDs in a special order with string  elements
    x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
    vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
    tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13] 
    reactor_model = params[4]; outflow_rate=reactor_model[0]; v_batch=reactor_model[1]; n_in=reactor_model[2]; SA=reactor_model[3]; v_R0=v_batch
    fat = params[5]; fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
    bat = params[6]; bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
    gat = params[7]; gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]
    dydt = np.zeros(9)

    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    Co = float(p_o/R/298.15) # in mol/m^3
    
    # reaction mechanism steps
    """
    R + 3X <=> RX + 3X
    RX + 2X <=> UX + HX + X
    UX + HX + X <=> VX + 2HX
    VX + 2HX <=> PX + 2HX
    PX + 2HX <=> P + 2HX + X
    P + 2HX + X <=> P + H2 + 3X
    or
    R + X <=> RX
    RX + X <=> UX + HX
    UX + X <=> VX + HX
    VX <=> PX
    PX <=> P + X
    2HX <=> H2 + 2X
    """

    # coverage and concentration definitions
    X  = y[0]    
    RX = y[1]
    UX = y[2]
    HX = y[3]
    PX = y[4]
    VX = y[5]
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
    rb4 = kb4 * n_P * X
    rf5 = kf5 * HX**2   
    rb5 = kb5 * n_H2 * (X**2)    
    """
    # case 1 when gas specie is in n/V, mol/m3.
    rf1 = kf1 * n_R/v_R0 * X
    rb1 = kb1 * RX  
    rf2 = kf2 * RX * X  
    rb2 = kb2 * UX * HX  
    rf3 = kf3 * UX * X   
    rb3 = kb3 * PX * HX  
    rf4 = kf4 * PX  
    rb4 = kb4 * n_P/v_R0 * X
    rf5 = kf5 * HX**2   
    rb5 = kb5 * n_H2/v_R0 * (X**2)
    rf6 = kf6 * VX  
    rb6 = kb6 * PX

    # ODEs for rate of change in surface specie coverages
    dydt[0] = 1*(-rf1 + rb1 - rf2 + rb2 - rf3 + rb3 + rf4 - rb4  + 2*rf5 - 2*rb5) # X in mol/m2/sec
    dydt[1] = 1*(rf1 - rb1 - rf2 + rb2) # RX in mol/m2/sec
    dydt[2] = 1*(rf2 - rb2 - rf3 + rb3) # UX in mol/m2/sec
    dydt[3] = 1*(rf2 - rb2 + rf3 - rb3 - 2*rf5 + 2*rb5) # HX in mol/m2/sec
    dydt[4] = 1*(rf3 - rb3 - rf4 + rb4 + rf6 - rb6) # PX in mol/m2/sec
    dydt[5] = 1*(-rf6 + rb6) # VX in mol/m2/sec
    #X = 1 - (RX+UX+HX+PX)semilogx

    # ODEs for rate of change in gas-phase concentration
    dydt[6] = SA*1e-6*(-rf1 + rb1) # R in mol/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2)  # H2 in mol/sec
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
    print('kf4=',kf4)    
    print('kb4=',kb4)
    print(kf4-kb4>0,kf4-kb4) 
    print('kf5=',kf5)    
    print('kb5=',kb5)
    print(kf5-kb5<0,kf5-kb5)

    # display of concentration(s) and coverage(s) for the intermediate & gas phase species 
    print('n_R=',n_R,'n_H2=',n_H2, 'n_P=',n_P, 'n_sum=',n_R+n_H2+n_P, 'n_sum11=',n_R+0.5*n_H2+0.5*n_P)
    print('RX=',RX, 'UX=',UX, 'HX=',HX, 'PX=',PX,'X=',1-RX-HX-UX-PX,'Xsum=',RX+HX+UX+PX)
    print()
    
    return dydt



































#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def solve_odes4444(n_R0,v_R0,mc0,C_R0,outflow_rate,T,tim,dydt,Gm,intID,specieID_list, S_x, fat, bat, gat):
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
  
    #Catalyst properties
    dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
    mc = mc0/1000 # mass of catalyst in g
    SA = dc * mc # surface area of catalyst in m^2
    sat_cov=18e-6 # mol/m2
    stand_cov=1.94e-7 # mol/m2
    sat_ML=sat_cov/sat_cov # in ML
    stand_ML=stand_cov/sat_cov # in ML
    #S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    #S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    #Input specifications for reactor model
    n_in = 0 # mol, amount of H2 charged in at t=0 
    #v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
    v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
    v_batch = v_initial # m3, volume the reactor/batch
    reactor_model = [outflow_rate, v_batch, n_in, SA]
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    p_o=1.01325e5 # standard pressure in Pa
    kB=const.Boltzmann # J/K
    N=const.Avogadro # mol^-1
    h=const.Planck # J.sam nie-7
    R=N*kB # const.gas_constant in J/mol/K

    energylist = Gm     # energy list
    int_id = intID;     # list of specie IDs  with string  elements
    rxn_specie_id = specieID_list     # list of selected specie IDs in a special order with string  elements
    x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
    vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
    tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13] 

    S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
    S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2   
    Co = float(p_o/R/298.15) # in mol/m^3

    fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
    bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
    gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]


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
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)*fat1/bat1
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T)*fat1 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)*fat2 
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)*bat2

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)*fat3 
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)*bat3

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T) 
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)*fat4 
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)*bat4  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)*bat5/fat5 
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T)*bat5  
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)*bat6/fat6 
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T)*bat6 
    kf5=kb5/keqq5*Co 

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    kf = [kf1,kf2,kf3,kf6,kf4,kf5]; kb = [kb1,kb2,kb3,kb6,kb4,kb5]
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    # Initial conditions
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim # total integration time in sec
    
    # construct ODE solver
    r = ode(dydt0).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=100, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, reactor_model, fat, bat, gat, kf, kb])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), 2)#
    yy = []
    tt = []
    time= []
    
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))        
        current=r.integrate(tnew)
        print('the reaction time (t_last) was ',tnew/60/60,' hours, with conversion of R as ',current[-3]/n_R0,' mol' )
        
        """
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))        
            current=r.integrate(tnew)
            print('the reaction time (t_max) was ',tnew/60/60,' hours')
            if current[-3]/n_R0<0.6:
                time.append(tnew/60/60)
                break
        print('the reaction time (t_max) was ',tnew/60/60,' hours')
        return tt, np.matrix(yy)     
        
        for x in xx:
            tnew = 10.0**x
            tt.append(tnew)
            yy.append(r.integrate(tnew))
        return tt, np.matrix(yy)
        """
        
    return tt, np.matrix(yy), kf, kb     

#=====================================================================================================
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS 
#=====================================================================================================
#=====================================================================================================


def dydt0(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]

    reactor_model = params[1]; outflow_rate=reactor_model[0]; v_batch=reactor_model[1]; n_in=reactor_model[2]; SA=reactor_model[3]; v_R0=v_batch
    fat = params[2]; fat1=fat[0]; fat2=fat[1]; fat3=fat[2]; fat4=fat[3]; fat5=fat[4]; fat6=fat[5]
    bat = params[3]; bat1=bat[0]; bat2=bat[1]; bat3=bat[2]; bat4=bat[3]; bat5=bat[4]; bat6=bat[5]
    gat = params[4]; gat1=gat[0]; gat2=gat[1]; gat3=gat[2]; gat4=gat[3]; gat5=gat[4]; gat6=gat[5]; gat7=gat[6]

    
    kf=params[5]
    kf1=kf[0]
    kf2=kf[1]
    kf3=kf[2]
    kf4=kf[3]
    kf5=kf[4]
    kf6=kf[5]
    kb=params[6]
    kb1=kb[0]
    kb2=kb[1]
    kb3=kb[2]
    kb4=kb[3]
    kb5=kb[4]
    kb6=kb[5]
    print()    

    dydt = np.zeros(9)

    
    # reaction mechanism steps
    """
    R + 3X <=> RX + 3X
    RX + 2X <=> UX + HX + X
    UX + HX + X <=> VX + 2HX
    VX + 2HX <=> PX + 2HX
    PX + 2HX <=> P + 2HX + X
    P + 2HX + X <=> P + H2 + 3X
    or
    R + X <=> RX
    RX + X <=> UX + HX
    UX + X <=> VX + HX
    VX <=> PX
    PX <=> P + X
    2HX <=> H2 + 2X
    """

    # coverage and concentration definitions
    X  = y[0]    
    RX = y[1]
    UX = y[2]
    HX = y[3]
    PX = y[4]
    VX = y[5]
    n_R = y[6]
    n_H2 = y[7]
    n_P = y[8]


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
    rb4 = kb4 * n_P * X
    rf5 = kf5 * HX**2   
    rb5 = kb5 * n_H2 * (X**2)    
    """
    # case 1 when gas specie is in n/V, mol/m3.
    rf1 = kf1 * n_R/v_R0 * X
    rb1 = kb1 * RX  
    rf2 = kf2 * RX * X  
    rb2 = kb2 * UX * HX  
    rf3 = kf3 * UX * X   
    rb3 = kb3 * PX * HX  
    rf4 = kf4 * PX  
    rb4 = kb4 * n_P/v_R0 * X
    rf5 = kf5 * HX**2   
    rb5 = kb5 * n_H2/v_R0 * (X**2)
    rf6 = kf6 * VX  
    rb6 = kb6 * PX

    # ODEs for rate of change in surface specie coverages
    dydt[0] = 1*(-rf1 + rb1 - rf2 + rb2 - rf3 + rb3 + rf4 - rb4  + 2*rf5 - 2*rb5) # X in mol/m2/sec
    dydt[1] = 1*(rf1 - rb1 - rf2 + rb2) # RX in mol/m2/sec
    dydt[2] = 1*(rf2 - rb2 - rf3 + rb3) # UX in mol/m2/sec
    dydt[3] = 1*(rf2 - rb2 + rf3 - rb3 - 2*rf5 + 2*rb5) # HX in mol/m2/sec
    dydt[4] = 1*(rf3 - rb3 - rf4 + rb4 + rf6 - rb6) # PX in mol/m2/sec
    dydt[5] = 1*(-rf6 + rb6) # VX in mol/m2/sec
    #X = 1 - (RX+UX+HX+PX)semilogx

    # ODEs for rate of change in gas-phase concentration
    dydt[6] = SA*1e-6*(-rf1 + rb1) # R in mol/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2)  # H2 in mol/sec
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
    print('kf4=',kf4)    
    print('kb4=',kb4)
    print(kf4-kb4>0,kf4-kb4) 
    print('kf5=',kf5)    
    print('kb5=',kb5)
    print(kf5-kb5<0,kf5-kb5)

    # display of concentration(s) and coverage(s) for the intermediate & gas phase species 
    print('n_R=',n_R,'n_H2=',n_H2, 'n_P=',n_P, 'n_sum=',n_R+n_H2+n_P, 'n_sum11=',n_R+0.5*n_H2+0.5*n_P)
    print('RX=',RX, 'UX=',UX, 'HX=',HX, 'PX=',PX,'X=',1-RX-HX-UX-PX,'Xsum=',RX+HX+UX+PX)
    print()

    return dydt


#=====================================================================================================
# CALLING FOR SOLUTION TO BE RUN
#=====================================================================================================
#=====================================================================================================

title1= '1-ol-oxacyclopentanol'
title2='propanol'
title3= 'cyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 


"""
TT = np.linspace(300, 1500, 10)  

def temp_effect1():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        xx,y = solve_odes(T,0.25,dydt1)      # <============= 
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1]))   
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(1)       # <=========
    plt.title(title1)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(11)       # <=========
    plt.title(title1)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx    

def temp_effect2():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        xx,y = solve_odes(T,0.25,dydt2)      # <============= 
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1]))   
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(2)       # <=========
    plt.title(title2)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(12)       # <=========
    plt.title(title2)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx   

def temp_effect3():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        xx,y = solve_odes(T,0.25,dydt3)      # <============= 
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1]))   
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(3)       # <=========
    plt.title(title3)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(13)       # <=========
    plt.title(title3)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx   

def temp_effect4():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        xx,y = solve_odes(T,0.25,dydt4)      # <============= 
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1]))   
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(4)       # <=========
    plt.title(title4)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(14)       # <=========
    plt.title(title4)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx   

def temp_effect5():      # <=============      
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    tlens=len(TT)
    i=0
    for i in range(tlens):
        T=TT[i]
        xx,y = solve_odes(T,0.25,dydt5)      # <============= 
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr.append(float(r[-1]))
        hh.append(float(h[-1]))
        pp.append(float(p[-1]))   
        xxcat.append(float(xcat[-1]/S_x)) 
        rrx.append(float(rx[-1]/S_x))
        uux.append(float(ux[-1]/S_x))
        hhx.append(float(hx[-1]/S_x))   
        vvx.append(float(vx[-1]/S_x))
        ppx.append(float(px[-1]/S_x))  
    max_p = max(pp)  # Find the maximum y value
    max_T = TT[pp.index(max_p)]  # Find the x value corresponding to the maximum y value
    print (max_T, round(max_p,5))
    print(rr,hh,pp,TT)
    plt.figure(5)       # <=========
    plt.title(title5)        # <=========    
    plt.plot(TT,rr, label='$n_R(T)$')
    plt.plot(TT,pp, label='$n_P(T)$')
    plt.plot(TT,hh, label='$n_{H2(T)}$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(15)       # <=========
    plt.title(title5)       # <=========
    plt.plot(TT,xxcat, label='$X$')
    plt.plot(TT,rrx, label='$RX$')
    plt.plot(TT,uux, label='$UX$')
    plt.plot(TT,hhx, label='$HX$')
    plt.plot(TT,vvx, label='$VX$')
    plt.plot(TT,ppx, label='$PX$')
    plt.xlabel('Temperature in K')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx   

rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=temp_effect1()
rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=temp_effect2()
rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=temp_effect3()
rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=temp_effect4()
rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=temp_effect5()

max_p1 = max(pp1)  # Find the maximum y value
max_p2 = max(pp2)  # Find the maximum y value
max_p3 = max(pp3)  # Find the maximum y value
max_p4 = max(pp4)  # Find the maximum y value
max_p5 = max(pp5)  # Find the maximum y value

max_T1 = TT[pp1.index(max_p1)]  # Find the x value corresponding to the maximum y value
max_T2 = TT[pp2.index(max_p2)]  # Find the x value corresponding to the maximum y value
max_T3 = TT[pp3.index(max_p3)]  # Find the x value corresponding to the maximum y value
max_T4 = TT[pp4.index(max_p4)]  # Find the x value corresponding to the maximum y value
max_T5 = TT[pp5.index(max_p5)]  # Find the x value corresponding to the maximum y value

print (title1,'Tmax1=', max_T1, 'Pmax1=', round(max_p1,5))
print (title2,'Tmax2=', max_T2, 'Pmax2=', round(max_p2,5))
print (title3,'Tmax3=', max_T3, 'Pmax3=', round(max_p3,5))
print (title4,'Tmax4=', max_T4, 'Pmax4=', round(max_p4,5))
print (title5,'Tmax5=', max_T5, 'Pmax5=', round(max_p5,5))

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

"""