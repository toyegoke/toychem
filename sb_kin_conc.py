from gn_ppt import G_m, int_id, S_m, H_m
from gn_thermop import specie_index
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const

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
# RELEVANT KINETIC INPUT DATA
#=====================================================================================================
#=====================================================================================================

#Reactor condition
T=145+273.15
p=1.01325e5 # standard pressure in Pa
kB=const.Boltzmann # J/K
N=const.Avogadro # mol^-1
h=const.Planck # J.sam ni
R=N*kB # const.gas_constant in J/mol/K
    
#Catalyst properties
dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
mc = 25/1000 # mass of catalyst in g
SA = dc * mc # surface area of catalyst in m^2
Asite = 1.96*1e-19 # site area in m2 (using the atomic radius 125e-12m or from reference: Tej & greeley, partial oxid. of meoh on moo3 010 dft)
sat_cov=18e-6 # mol/m2
stand_cov=1.94e-7 # mol/m2
sat_ML=sat_cov/sat_cov # in ML
stand_ML=stand_cov/sat_cov # in ML
S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
S_o_in = (1.94*1e-7)# mol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2
S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
S_o = S_o_in/(1e-6) # umol/m^2 # np.exp(1/3)*(Co)**(2/3) # mol/m^2

#Feed parameters
Co = p/R/298.15 # in mol/m^3
C_R0 = 950 # in mol/m3 (0.95 M)
v_R0 = 45*1e-6 # in m3 (45 mL)
n_R0 = v_R0*C_R0 # in mol

#Input specifications for reactor model
n_in = 0 # mol, amount of H2 charged in at t=0
outflow_rate = 5e-7 # m3/sec, volumetric flowrate for the H2 withdrawal 
#v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
v_batch = v_initial # m3, volume the reactor/batch

#=====================================================================================================
#=====================================================================================================
# ODE MODEL SOLVER
#=====================================================================================================
#=====================================================================================================

def solve_odes(T,tim,dydt,xC):
    """
    Time-integrate chemo-kinetic system
    """
    # initial conditions
    n_R00=n_R0*(1+xC) # in mol
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R00, 0, 0]
    t0 = 0
    t1 = tim*60*60 # total integration time

    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=500, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    xx = np.linspace(0, np.log10(t1), int((np.log10(t1) + 12.0)*10))
    yy = []
    tt = []
    time= []
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

    """
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))
    return tt, np.matrix(yy)
    """

#=====================================================================================================
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR 2-PROPANOL
#=====================================================================================================
#=====================================================================================================

def dydt1(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
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

    #Free Gibb (G) energy for the species
    Xg=G_m[specie_index('X',int_id)]
    Rg=G_m[specie_index('Ra',int_id)]
    RXg=G_m[specie_index('RaX',int_id)]
    UXg=G_m[specie_index('UaX',int_id)]
    HXg=G_m[specie_index('HX',int_id)]
    VXg=G_m[specie_index('VaX',int_id)]
    PXg=G_m[specie_index('PaX',int_id)]
    H2g=G_m[specie_index('H2',int_id)]
    Pg=G_m[specie_index('Pa',int_id)]
    TS1g=G_m[specie_index('TS1a',int_id)]
    TS2g=G_m[specie_index('TS2a',int_id)]
    TSRg=G_m[specie_index('TSRa',int_id)]
    TSPg=G_m[specie_index('TSPa',int_id)]
    TSH2g=G_m[specie_index('TSH2',int_id)]

    # calculate all reaction rate constants
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T) 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T)
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T) 
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T) 
    kf5=kb5/keqq5*Co

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
    dydt[6] = SA*1e-6*(-rf1 + rb1) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) #+ outflow_rate/v_batch*(Cin - y[7]) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec

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
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR CYCLOPENTANOL
#=====================================================================================================
#=====================================================================================================

def dydt2(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
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

    #Free Gibb (G) energy for the species
    Xg=G_m[specie_index('X',int_id)]
    Rg=G_m[specie_index('Rb',int_id)]
    RXg=G_m[specie_index('RbX',int_id)]
    UXg=G_m[specie_index('UbX',int_id)]
    HXg=G_m[specie_index('HX',int_id)]
    VXg=G_m[specie_index('VbX',int_id)]
    PXg=G_m[specie_index('PbX',int_id)]
    H2g=G_m[specie_index('H2',int_id)]
    Pg=G_m[specie_index('Pb',int_id)]
    TS1g=G_m[specie_index('TS1b',int_id)]
    TS2g=G_m[specie_index('TS2b',int_id)]
    TSRg=G_m[specie_index('TSRb',int_id)]
    TSPg=G_m[specie_index('TSPb',int_id)]
    TSH2g=G_m[specie_index('TSH2',int_id)]

    # calculate all reaction rate constants
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T) 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T)
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T) 
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T) 
    kf5=kb5/keqq5*Co

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
    dydt[6] = SA*1e-6*(-rf1 + rb1) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) #+ outflow_rate/v_batch*(Cin - y[7]) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec

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
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR 1-OXACYCLOPENTANOL
#=====================================================================================================
#=====================================================================================================

def dydt3(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
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

    #Free Gibb (G) energy for the species
    Xg=G_m[specie_index('X',int_id)]
    Rg=G_m[specie_index('R',int_id)]
    RXg=G_m[specie_index('RX',int_id)]
    UXg=G_m[specie_index('UX',int_id)]
    HXg=G_m[specie_index('HX',int_id)]
    VXg=G_m[specie_index('VX',int_id)]
    PXg=G_m[specie_index('PX',int_id)]
    H2g=G_m[specie_index('H2',int_id)]
    Pg=G_m[specie_index('P',int_id)]
    TS1g=G_m[specie_index('TS1',int_id)]
    TS2g=G_m[specie_index('TS2',int_id)]
    TSRg=G_m[specie_index('TSR',int_id)]
    TSPg=G_m[specie_index('TSP',int_id)]
    TSH2g=G_m[specie_index('TSH2',int_id)]

    # calculate all reaction rate constants
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T) 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T)
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T) 
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T) 
    kf5=kb5/keqq5*Co

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
    dydt[6] = SA*1e-6*(-rf1 + rb1) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) #+ outflow_rate/v_batch*(Cin - y[7]) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec

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
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR 2-OXACYCLOPENTANOL
#=====================================================================================================
#=====================================================================================================

def dydt4(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
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

    #Free Gibb (G) energy for the species
    Xg=G_m[specie_index('X',int_id)]
    Rg=G_m[specie_index('Rc',int_id)]
    RXg=G_m[specie_index('RcX',int_id)]
    UXg=G_m[specie_index('UcX',int_id)]
    HXg=G_m[specie_index('HX',int_id)]
    VXg=G_m[specie_index('VcX',int_id)]
    PXg=G_m[specie_index('PcX',int_id)]
    H2g=G_m[specie_index('H2',int_id)]
    Pg=G_m[specie_index('Pc',int_id)]
    TS1g=G_m[specie_index('TS1c',int_id)]
    TS2g=G_m[specie_index('TS2c',int_id)]
    TSRg=G_m[specie_index('TSRc',int_id)]
    TSPg=G_m[specie_index('TSPc',int_id)]
    TSH2g=G_m[specie_index('TSH2',int_id)]

    # calculate all reaction rate constants
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T) 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T)
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T) 
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T)
    kf5=kb5/keqq5*Co

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
    dydt[6] = SA*1e-6*(-rf1 + rb1) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) #+ outflow_rate/v_batch*(Cin - y[7]) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec

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
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR 2-ADAMANTANOL (FFCC)
#=====================================================================================================
#=====================================================================================================

def dydt5(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
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

    #Free Gibb (G) energy for the species
    Xg=G_m[specie_index('X',int_id)]
    Rg=G_m[specie_index('Rd',int_id)]
    RXg=G_m[specie_index('RdX',int_id)]
    UXg=G_m[specie_index('UdX',int_id)]
    HXg=G_m[specie_index('HX',int_id)]
    VXg=G_m[specie_index('VdX',int_id)]
    PXg=G_m[specie_index('PdX',int_id)]
    H2g=G_m[specie_index('H2',int_id)]
    Pg=G_m[specie_index('Pd',int_id)]
    TS1g=G_m[specie_index('TS1d',int_id)]
    TS2g=G_m[specie_index('TS2d',int_id)]
    TSRg=G_m[specie_index('TSRd',int_id)]
    TSPg=G_m[specie_index('TSPd',int_id)]
    TSH2g=G_m[specie_index('TSH2',int_id)]

    # calculate all reaction rate constants
    print('R+X==>RX')
    keqq1=kineticparameter.keq(Xg+Rg,RXg,T)
    kf1=kineticparameter.ksu(TSRg,Rg+Xg,S_o,1,1,T) 
    kb1=kf1/keqq1*Co

    print('RX+X==>UX+HX')
    keqq2=kineticparameter.keq(RXg+Xg,UXg+HXg,T)   
    kf2=kineticparameter.ksu(TS1g+Xg,RXg+Xg,S_o,0,2,T)
    kb2=kineticparameter.ksu(TS1g+Xg,UXg+HXg,S_o,0,2,T)

    print('UX+X==>VX+HX')
    keqq3=kineticparameter.keq(UXg+Xg,VXg+HXg,T)
    kf3=kineticparameter.ksu(TS2g+Xg,UXg+Xg,S_o,0,2,T)
    kb3=kineticparameter.ksu(TS2g+Xg,VXg+HXg,S_o,0,2,T)

    print('VX==>PX')
    keqq3=kineticparameter.keq(PXg,VXg,T)
    kf6=kineticparameter.ksu(PXg,VXg,S_o,0,1,T)
    kb6=kineticparameter.ksu(PXg,PXg,S_o,0,1,T)  

    print('PX==>P+X')
    keqq4=kineticparameter.keq(Pg+Xg,PXg,T)
    kb4=kineticparameter.ksu(TSPg,Pg+Xg,S_o,1,1,T) 
    kf4=kb4/keqq4*Co
    
    print('2HX==>H2+2X')
    keqq5=kineticparameter.keq(H2g+2*Xg,2*HXg,T)
    kb5=kineticparameter.ksu(TSH2g,H2g+2*Xg,S_o,1,1,T)  
    kf5=kb5/keqq5*Co

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
    dydt[6] = SA*1e-6*(-rf1 + rb1) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) #+ outflow_rate/v_batch*(Cin - y[7]) #1/N/nc#SA*sat_surf_conc/vg#1/N/vg/sat_surf_conc #SA/vg #nc/vg # R in mol/m^3/sec or M/sec

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
#=====================================================================================================
# CALLING FOR SOLUTION TO BE RUN
#=====================================================================================================
#=====================================================================================================

title2='propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC'
xCC = np.linspace(-0.95, 0.95, 10)  # fraction use to study the change in nR in decimal
nR = n_R0*(1+xCC) # the actual changes in nR evaluated in moles


def conc_effect1():   # <=========
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    conclens=len(xCC)
    i=0
    for i in range(conclens):
        T=800 # in K
        xC=xCC[i]
        xx,y = solve_odes(T,0.25,dydt1,xC) # <=========
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr0=n_R0*(1+xC)
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
    max_nR = n_R0*(1+xCC[pp.index(max_p)])  # Find the x value corresponding to the maximum y value
    print (max_nR, round(max_p,5))
    print(rr,hh,pp,xCC)
    plt.figure(1)       # <=========
    plt.title(title1)        # <=========
    plt.plot(nR,rr, label='$R(unconv)$')
    plt.plot(nR,pp, label='$P(prod)$')
    plt.plot(nR,hh, label='$H(prod)$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(11)       # <=========
    plt.title(title1)       # <=========
    plt.plot(nR,xxcat, label='$X$')
    plt.plot(nR,rrx, label='$RX$')
    plt.plot(nR,uux, label='$UX$')
    plt.plot(nR,hhx, label='$HX$')
    plt.plot(nR,vvx, label='$VX$')
    plt.plot(nR,ppx, label='$PX$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx


def conc_effect2():   # <=========
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    conclens=len(xCC)
    i=0
    for i in range(conclens):
        T=800 # in K
        xC=xCC[i]
        xx,y = solve_odes(T,0.25,dydt2,xC) # <=========
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr0=n_R0*(1+xC)
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
    max_nR = n_R0*(1+xCC[pp.index(max_p)])  # Find the x value corresponding to the maximum y value
    print (max_nR, round(max_p,5))
    print(rr,hh,pp,xCC)
    plt.figure(2)       # <=========
    plt.title(title2)        # <=========
    plt.plot(nR,rr, label='$R(unconv)$')
    plt.plot(nR,pp, label='$P(prod)$')
    plt.plot(nR,hh, label='$H(prod)$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(12)       # <=========
    plt.title(title2)       # <=========
    plt.plot(nR,xxcat, label='$X$')
    plt.plot(nR,rrx, label='$RX$')
    plt.plot(nR,uux, label='$UX$')
    plt.plot(nR,hhx, label='$HX$')
    plt.plot(nR,vvx, label='$VX$')
    plt.plot(nR,ppx, label='$PX$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx

def conc_effect3():   # <=========
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    conclens=len(xCC)
    i=0
    for i in range(conclens):
        T=800 # in K
        xC=xCC[i]
        xx,y = solve_odes(T,0.25,dydt3,xC) # <=========
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr0=n_R0*(1+xC)
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
    max_nR = n_R0*(1+xCC[pp.index(max_p)])  # Find the x value corresponding to the maximum y value
    print (max_nR, round(max_p,5))
    print(rr,hh,pp,xCC)
    plt.figure(3)       # <=========
    plt.title(title3)        # <=========
    plt.plot(nR,rr, label='$R(unconv)$')
    plt.plot(nR,pp, label='$P(prod)$')
    plt.plot(nR,hh, label='$H(prod)$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(13)       # <=========
    plt.title(title3)       # <=========
    plt.plot(nR,xxcat, label='$X$')
    plt.plot(nR,rrx, label='$RX$')
    plt.plot(nR,uux, label='$UX$')
    plt.plot(nR,hhx, label='$HX$')
    plt.plot(nR,vvx, label='$VX$')
    plt.plot(nR,ppx, label='$PX$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx

def conc_effect4():   # <=========
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    conclens=len(xCC)
    i=0
    for i in range(conclens):
        T=800 # in K
        xC=xCC[i]
        xx,y = solve_odes(T,0.25,dydt4,xC) # <=========
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr0=n_R0*(1+xC)
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
    max_nR = n_R0*(1+xCC[pp.index(max_p)])  # Find the x value corresponding to the maximum y value
    print (max_nR, round(max_p,5))
    print(rr,hh,pp,xCC)
    plt.figure(4)       # <=========
    plt.title(title4)        # <=========
    plt.plot(nR,rr, label='$R(unconv)$')
    plt.plot(nR,pp, label='$P(prod)$')
    plt.plot(nR,hh, label='$H(prod)$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(14)       # <=========
    plt.title(title4)       # <=========
    plt.plot(nR,xxcat, label='$X$')
    plt.plot(nR,rrx, label='$RX$')
    plt.plot(nR,uux, label='$UX$')
    plt.plot(nR,hhx, label='$HX$')
    plt.plot(nR,vvx, label='$VX$')
    plt.plot(nR,ppx, label='$PX$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx

def conc_effect5():   # <=========
    rr = []
    hh = []
    pp = []
    xxcat = []
    rrx = []
    uux = []
    hhx = []
    vvx = []
    ppx = []
    conclens=len(xCC)
    i=0
    for i in range(conclens):
        T=800 # in K
        xC=xCC[i]
        xx,y = solve_odes(T,0.25,dydt5,xC) # <=========
        xcat=y[:,0]#*1e-6
        rx=y[:,1]#*1e-6
        ux=y[:,2]#*1e-6
        hx=y[:,3]#*1e-6
        vx=y[:,4]#*1e-6
        px=y[:,5]#*1e-6
        r=y[:,6]#*1e-6
        h=y[:,7]#*1e-6
        p=y[:,8]#*1e-6
        rr0=n_R0*(1+xC)
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
    max_nR = n_R0*(1+xCC[pp.index(max_p)])  # Find the x value corresponding to the maximum y value
    print (max_nR, round(max_p,5))
    print(rr,hh,pp,xCC)
    plt.figure(5)       # <=========
    plt.title(title5)        # <=========
    plt.plot(nR,rr, label='$R(unconv)$')
    plt.plot(nR,pp, label='$P(prod)$')
    plt.plot(nR,hh, label='$H(prod)$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Amount of Substance in mol')
    plt.legend(loc=2,prop={'size':8})
    plt.figure(15)       # <=========
    plt.title(title5)       # <=========
    plt.plot(nR,xxcat, label='$X$')
    plt.plot(nR,rrx, label='$RX$')
    plt.plot(nR,uux, label='$UX$')
    plt.plot(nR,hhx, label='$HX$')
    plt.plot(nR,vvx, label='$VX$')
    plt.plot(nR,ppx, label='$PX$')
    plt.xlabel('Change in Feed Qty, R in mol')
    plt.ylabel('Coverage Fraction')
    plt.legend(loc=5,prop={'size':8})
    return rr,hh,pp,xxcat,rrx,uux,hhx,vvx,ppx


rr1,hh1,pp1,x1,rx1,ux1,hx1,vx1,px1=conc_effect1()
rr2,hh2,pp2,x2,rx2,ux2,hx2,vx2,px2=conc_effect2()
rr3,hh3,pp3,x3,rx3,ux3,hx3,vx3,px3=conc_effect3()
rr4,hh4,pp4,x4,rx4,ux4,hx4,vx4,px4=conc_effect4()
rr5,hh5,pp5,x5,rx5,ux5,hx5,vx5,px5=conc_effect5()


 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx