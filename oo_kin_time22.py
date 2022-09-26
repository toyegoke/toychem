from ppt import G_m, int_id, S_m, H_m
from thermop import specie_index
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
        K=1 # transmission coefficient factor
        p=1.01325e5 # standard pressure in Pa
        kB=const.Boltzmann
        N=const.Avogadro
        h=const.Planck
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
T=550 # 145+273.15
p=1.01325e5 # standard pressure in Pa
kB=const.Boltzmann # J/K
N=const.Avogadro # mol^-1
h=const.Planck # J.sam ni
R=N*kB # const.gas_constant in J/mol/K
    
#Catalyst properties
dc = 29#(29+26+39+39+27+36+27+35+47+39+28+22+28+23)/14 # specific surface area of catalyst in m^2/g
mc = 40/1000 # 25/1000 # mass of catalyst in g
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
C_R0 = 998   # 950 # in mol/m3 (0.95 M)
v_R0 = 45*1e-6 # in m3 (45 mL)
n_R0 = v_R0*C_R0 # in mol

#Input specifications for reactor model
n_in = 0 # mol, amount of H2 charged in at t=0
outflow_rate = 14.27e-7  #5e-7 # m3/sec, volumetric flowrate for the H2 withdrawal 
#v_initial = n_R0*R*T/p # m3, volume the reactor/batch when it is dependent on temp (T)
v_initial = 2*v_R0  # m3, volume the reactor/batch when it is dependent on temp (T)
v_batch = v_R0 # v_initial # m3, volume the reactor/batch

#=====================================================================================================
#=====================================================================================================
# ODE MODEL SOLVER
#=====================================================================================================
#=====================================================================================================

def solve_odes(T,tim,dydt):
    """
    Time-integrate chemo-kinetic system
    """
    # initial conditions
    n_R0 # in mol
    S_x0 = S_x
    y0 = [S_x0, 0, 0, 0, 0, 0, n_R0, 0, 0]
    t0 = 0
    t1 = tim*3600 # total integration time

    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf', 
           atol=1e-8, rtol=1e-8, nsteps=500, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T])

    # integrate on a logaritmic scale
    # xx = np.linspace(0, 1000, 100)
    #xx = np.linspace(0, np.log10(t1), int((np.log10(t1) + 12.0)*10))
    xx = np.linspace(0, np.log10(t1), 50)
    #xx = np.linspace(0, t1, 20)
    yy = []
    tt = []
    time= []
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))
    return tt, np.matrix(yy)


#=====================================================================================================
#=====================================================================================================
# MICRO-KINECTIC & REACTOR MODELS FOR 2-PROPANOL
#=====================================================================================================
#=====================================================================================================

def dydt1(t, y, params):
    """
    Set of ordinary differential equations
    """
    kB=const.Boltzmann
    N=const.Avogadro
    h=const.Planck
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
    kB=const.Boltzmann
    N=const.Avogadro
    h=const.Planck
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
    kB=const.Boltzmann
    N=const.Avogadro
    h=const.Planck
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
    kB=const.Boltzmann
    N=const.Avogadro
    h=const.Planck
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
    kB=const.Boltzmann
    N=const.Avogadro
    h=const.Planck
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
    dydt[6] = SA*1e-6*(-rf1 + rb1) ## R in mol/m^3/sec or M/sec 
    dydt[7] = SA*1e-6*(rf5 - rb5) + outflow_rate/v_batch*(n_in - n_H2) ## H2 in mol/m^3/sec or M/sec
    dydt[8] = SA*1e-6*(rf4 - rb4) ## P in mol/m^3/sec or M/sec

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

def sol(moleculename,temperature,time,odeset):
    rr0=1 # n_R0 # initial amount of R 
    time_all = time # hr
    T1 = temperature
    dyd = odeset
    
    xx1,y1 = solve_odes(T1,time_all,dyd)
    xl1=len(xx1)
    x1=[]
    for i in range(xl1):
        x1.append(xx1[i]/3600) # sec to hr # that convert it from seconds into hours
    XX1=1-(y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5])
    rr1=y1[:,6]#*1e-6
    hh1=y1[:,7]#*1e-6
    pp1=y1[:,8]#*1e-6
    
    fig1=plt.figure()
    plt.plot(x1, rr1/rr0, label='$n_R$')
    plt.plot(x1, (pp1-hh1)/rr0, label='n_H2(out)')
    plt.plot(x1, (hh1)/rr0, label='n_H2(in)')
    plt.plot(x1, pp1/rr0, label='$n_P$')
    plt.xlabel('Time (t) in h')
    plt.ylabel('Amount (n) in mol')
    plt.title(moleculename)
    plt.legend()
    fig1.savefig('/home/toyegoke/ENS_kinetics/Final/'+moleculename+'_plot_data1_for_RDSk(P_G).png')
 
    fig6=plt.figure()
    plt.plot(x1, y1[:,0]/S_x, label='c_X')
    plt.plot(x1, y1[:,1]/S_x, label='c_RX')
    plt.plot(x1, y1[:,2]/S_x, label='c_UX')
    plt.plot(x1, y1[:,3]/S_x, label='C_HX')
    plt.plot(x1, y1[:,4]/S_x, label='c_PX')
    plt.plot(x1, y1[:,5]/S_x, label='c_VX')
    plt.xlabel('Time (t) in h')
    plt.ylabel('Coverage (Si/Sx) in ML')
    plt.legend()  
    plt.title(moleculename)
    fig6.savefig('/home/toyegoke/ENS_kinetics/Final/'+moleculename+'_plot6_data_for_RDSk(P_G).png')  
    
    print('2-propanol results')
    print('==================================')
    print('n(P)=',pp1[-1],'mol')
    print('n(H2)=',hh1[-1],'mol')
    print('n(R)=',rr1[-1],'mol')
    print('Y(P)=',pp1[-1]/rr0)
    print('Y(H2)=',hh1[-1]/rr0)
    print('X(R)=',(rr0-rr1[-1])/rr0)
    print('==================================')
    
    print('Overall Summary Results')
    print('==================================')
    print('X(R)=',(rr0-rr1[-1])/rr0,'t(h)=',x1[-1]/3600,moleculename)
    print('==================================')
    
    return
    

#=====================================================================================================
# CALLING FOR SOLUTION TO BE RUN
#=====================================================================================================
#=====================================================================================================

title1= '1-ol-oxacyclopentanol'
title2= 'propanol'
title3= 'cyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 

titl = [title1,title2,title3,title4,title5]
pt = [dydt1,dydt2,dydt3,dydt4,dydt5]
ranP = len(titl)
i=0

for i in range(ranP):
    ppt=pt[i]
    title=titl[i]
    sol(title,550,14.49/60,ppt)#14.49/60





"""
rr0=1 # n_R0 # initial amount of R 
time_all = 24 # hr
#T1=1500 # K
#T2=700 # K
#T3=567 # K
#T4=1500 # K
#T5=1233 # K

T1=550
T2=T1; T3=T1; T4=T1; T5=T1
 
# solution for 2-propanol
xx1,y1 = solve_odes(T1,time_all,dydt1)
xl1=len(xx1)
x1=[]
for i in range(xl1):
    x1.append(xx1[i]/3600) # sec to hr # that convert it from seconds into hours
XX1=1-(y1[:,1]+y1[:,2]+y1[:,3]+y1[:,4]+y1[:,5])
rr1=y1[:,6]#*1e-6
hh1=y1[:,7]#*1e-6
pp1=y1[:,8]#*1e-6

# solution for cyclopentanol
xx2,y2 = solve_odes(T2,time_all,dydt2)
xl2=len(xx2)
x2=[]
for i in range(xl2):
    x2.append(xx2[i]/3600) # sec to hr # that convert it from seconds into hours
XX2=1-(y2[:,1]+y2[:,2]+y2[:,3]+y2[:,4]+y2[:,5])
rr2=y2[:,6]#*1e-6
hh2=y2[:,7]#*1e-6
pp2=y2[:,8]#*1e-6

# solution for 1-oxacyclopentanol
xx3,y3 = solve_odes(T3,time_all,dydt3)
xl3=len(xx3)
x3=[]
for i in range(xl3):
    x3.append(xx3[i]/3600) # sec to hr # that convert it from seconds into hours
XX3=1-(y3[:,1]+y3[:,2]+y3[:,3]+y3[:,4]+y3[:,5])
rr3=y3[:,6]#*1e-6
hh3=y3[:,7]#*1e-6
pp3=y3[:,8]#*1e-6

# solution for 2-oxacyclopentanol
xx4,y4 = solve_odes(T4,time_all,dydt4)
xl4=len(xx4)
x4=[]
for i in range(xl4):
    x4.append(xx4[i]/3600) # sec to hr # that convert it from seconds into hours
XX4=1-(y4[:,1]+y4[:,2]+y4[:,3]+y4[:,4]+y4[:,5])
rr4=y4[:,6]#*1e-6
hh4=y4[:,7]#*1e-6
pp4=y4[:,8]#*1e-6

# solution for 2-adamantanol
xx5,y5 = solve_odes(T5,time_all,dydt5)
xl5=len(xx5)
x5=[]
for i in range(xl5):
    x5.append(xx5[i]/3600) # sec to hr # that convert it from seconds into hours
XX5=1-(y5[:,1]+y5[:,2]+y5[:,3]+y5[:,4]+y5[:,5])
rr5=y5[:,6]#*1e-6
hh5=y5[:,7]#*1e-6
pp5=y5[:,8]#*1e-6




# solution for 2-propanol
fig1=plt.figure(4)
plt.plot(x1, rr1/rr0, label='$n_R(2-propanol)$')
plt.plot(x1, hh1/rr0, label='n_H2(out)')
plt.plot(x1, pp1/rr0, label='$n_P(acetone)$')
plt.xlabel('Time (t) in h')
plt.ylabel('Amount (n) in mol')
plt.legend()
fig1.savefig('/home/toyegoke/ENS_kinetics/Final/_plot_data1_for_RDSk(P_G).png')

# solution for cyclopentanol
fig2=plt.figure(5)
plt.plot(x2, rr2/rr0, label='$n_R(cyclopentanol)$')
plt.plot(x2, hh2/rr0, label='n_H2(out)')
plt.plot(x2, pp2/rr0, label='$n_P(cyclopentanone)$')
plt.xlabel('Time (t) in h')
plt.ylabel('Amount (n) in mol')
plt.legend()
fig2.savefig('/home/toyegoke/ENS_kinetics/Final/_plot_data2_for_RDSk(P_G).png')

# solution for 1-oxacyclopentanol
fig3=plt.figure(6)
plt.plot(x3, rr3/rr0, label='$n_R(1-oxacyclopentanol)$')
plt.plot(x3, hh3/rr0, label='n_H2(out)')
plt.plot(x3, pp3/rr0, label='$n_P(1-oxacyclopentanone)$')
plt.xlabel('Time (t) in h')
plt.ylabel('Amount (n) in mol')
plt.legend()
fig3.savefig('/home/toyegoke/ENS_kinetics/Final/_plot3_data_for_RDSk(P_G).png')

# solution for 2-oxacyclopentanol
fig4=plt.figure(7)
plt.plot(x4, rr4/rr0, label='$n_R(2-oxacyclopentanol)$')
plt.plot(x4, hh4/rr0, label='n_H2(out)')
plt.plot(x4, pp4/rr0, label='$n_P(2-oxacyclopentanone)$')
plt.xlabel('Time (t) in h')
plt.ylabel('Amount (n) in mol')
plt.legend()
fig4.savefig('/home/toyegoke/ENS_kinetics/Final/_plot4_data_for_RDSk(P_G).png')

# solution for 2-adamantanol
fig5=plt.figure(8)
plt.plot(x5, rr5/rr0, label='$n_R(2-adamantanol)$')
plt.plot(x5, hh5/rr0, label='n_H2(out)')
plt.plot(x5, pp5/rr0, label='$n_P(2-adamantanone)$')
plt.xlabel('Time (t) in h')
plt.ylabel('Amount (n) in mol')
plt.legend()
fig5.savefig('/home/toyegoke/ENS_kinetics/Final/_plot5_data_for_RDSk(P_G).png')



# solution for 2-propanol
fig6=plt.figure(10)
plt.plot(x1, y1[:,0]/S_x, label='c_X')
plt.plot(x1, y1[:,1]/S_x, label='c_RX')
plt.plot(x1, y1[:,2]/S_x, label='c_UX')
plt.plot(x1, y1[:,3]/S_x, label='C_HX')
plt.plot(x1, y1[:,4]/S_x, label='c_PX')
plt.plot(x1, y1[:,5]/S_x, label='c_VX')
plt.xlabel('Time (t) in h (2-propanol Case)')
plt.ylabel('Coverage (Si/Sx) in ML')
plt.legend()  
fig6.savefig('/home/toyegoke/ENS_kinetics/Final/_plot6_data_for_RDSk(P_G).png')  

# solution for cyclopentanol
fig7=plt.figure(11)
plt.plot(x2, y2[:,0]/S_x, label='c_X')
plt.plot(x2, y2[:,1]/S_x, label='c_RX')
plt.plot(x2, y2[:,2]/S_x, label='c_UX')
plt.plot(x2, y2[:,3]/S_x, label='C_HX')
plt.plot(x2, y2[:,4]/S_x, label='c_PX')
plt.plot(x2, y2[:,5]/S_x, label='c_VX')
plt.xlabel('Time (t) in h (cyclopentanol Case)')
plt.ylabel('Coverage (Si/Sx) in ML')
plt.legend()    
fig7.savefig('/home/toyegoke/ENS_kinetics/Final/_plot7_data_for_RDSk(P_G).png')

# solution for 1-oxacyclopentanol
fig8=plt.figure(12)
plt.plot(x3, y3[:,0]/S_x, label='c_X')
plt.plot(x3, y3[:,1]/S_x, label='c_RX')
plt.plot(x3, y3[:,2]/S_x, label='c_UX')
plt.plot(x3, y3[:,3]/S_x, label='C_HX')
plt.plot(x3, y3[:,4]/S_x, label='c_PX')
plt.plot(x3, y3[:,5]/S_x, label='c_VX')
plt.xlabel('Time (t) in h (1-oxacyclopentanol Case)')
plt.ylabel('Coverage (Si/Sx) in ML')
plt.legend()    
fig8.savefig('/home/toyegoke/ENS_kinetics/Final/_plot8_data_for_RDSk(P_G).png')

# solution for 2-oxacyclopentanol
fig9=plt.figure(13)
plt.plot(x4, y4[:,0]/S_x, label='c_X')
plt.plot(x4, y4[:,1]/S_x, label='c_RX')
plt.plot(x4, y4[:,2]/S_x, label='c_UX')
plt.plot(x4, y4[:,3]/S_x, label='C_HX')
plt.plot(x4, y4[:,4]/S_x, label='c_PX')
plt.plot(x4, y4[:,5]/S_x, label='c_VX')
plt.xlabel('Time (t) in h (2-oxacyclopentanol Case)')
plt.ylabel('Coverage (Si/Sx) in ML')
plt.legend()    
fig9.savefig('/home/toyegoke/ENS_kinetics/Final/_plot9_data_for_RDSk(P_G).png')

# solution for 2-adamantanol
fig10=plt.figure(14)
plt.plot(x5, y5[:,0]/S_x, label='c_X')
plt.plot(x5, y5[:,1]/S_x, label='c_RX')
plt.plot(x5, y5[:,2]/S_x, label='c_UX')
plt.plot(x5, y5[:,3]/S_x, label='C_HX')
plt.plot(x5, y5[:,4]/S_x, label='c_PX')
plt.plot(x5, y5[:,5]/S_x, label='c_VX')
plt.xlabel('Time (t) in h (2-adamantanol Case)')
plt.ylabel('Coverage (Si/Sx) in ML')
plt.legend()    
fig10.savefig('/home/toyegoke/ENS_kinetics/Final/_plot10_data_for_RDSk(P_G).png')

# solution for 2-propanol
fig11=plt.figure(20)

#lt.plot(x1, (pp1 -0) / (rr0), label='Y(P,acetone)')
plt.plot(x1, (rr0-rr1) / (rr0), label='X(R,2-propanol)')

#lt.plot(x2, (pp2 -0) / (rr0), label='Y(P,cyclopentanone)')
plt.plot(x2, (rr0-rr2) / (rr0), label='X(R,cyclopentanol)')

#lt.plot(x3, (pp3 -0) / (rr0), label='Y(P,1-oxacyclopentanone)')
plt.plot(x3, (rr0-rr3) / (rr0), label='X(R,1-oxacyclopentanol)')

#lt.plot(x4, (pp4 -0) / (rr0), label='Y(P,2-oxacyclopentanone)')
plt.plot(x4, (rr0-rr4) / (rr0), label='X(R,2-oxacyclopentanol)')

#lt.plot(x5, (pp5 -0) / (rr0), label='Y(P,2-adamantanone)')
plt.plot(x5, (rr0-rr5) / (rr0), label='X(R,2-adamantanol)')

plt.xlabel('Time (t) in h')
plt.ylabel('X(R)') #Y(P)/
plt.legend()
#plt.show()
fig11.savefig('/home/toyegoke/ENS_kinetics/Final/_plot_data11_for_RDSk(P_G).png')

print('2-propanol results')
print('==================================')
print('n(P)=',pp1[-1],'mol')
print('n(H2)=',hh1[-1],'mol')
print('n(R)=',rr1[-1],'mol')
print('Y(P)=',pp1[-1]/rr0)
print('Y(H2)=',hh1[-1]/rr0)
print('X(R)=',(rr0-rr1[-1])/rr0)
print('==================================')

print('cyclopentanol results')
print('==================================')
print('n(P)=',pp2[-1],'mol')
print('n(H2)=',hh2[-1],'mol')
print('n(R)=',rr2[-1],'mol')
print('Y(P)=',pp2[-1]/rr0)
print('Y(H2)=',hh2[-1]/rr0)
print('X(R)=',(rr0-rr2[-1])/rr0)
print('==================================')

print('1-oxacyclopentanol results')
print('==================================')
print('n(P)=',pp3[-1],'mol')
print('n(H2)=',hh3[-1],'mol')
print('n(R)=',rr3[-1],'mol')
print('Y(P)=',pp3[-1]/rr0)
print('Y(H2)=',hh3[-1]/rr0)
print('X(R)=',(rr0-rr3[-1])/rr0)
print('==================================')

print('2-oxacyclopentanol results')
print('==================================')
print('n(P)=',pp4[-1],'mol')
print('n(H2)=',hh4[-1],'mol')
print('n(R)=',rr4[-1],'mol')
print('Y(P)=',pp4[-1]/rr0)
print('Y(H2)=',hh4[-1]/rr0)
print('X(R)=',(rr0-rr4[-1])/rr0)
print('==================================')

print('2-adamantanol results')
print('==================================')
print('n(P)=',pp5[-1],'mol')
print('n(H2)=',hh5[-1],'mol')
print('n(R)=',rr5[-1],'mol')
print('Y(P)=',pp5[-1]/rr0)
print('Y(H2)=',hh5[-1]/rr0)
print('X(R)=',(rr0-rr5[-1])/rr0)
print('==================================')

print('Overall Summary Results')
print('==================================')
print('X(R)=',(rr0-rr1[-1])/rr0,'t(h)=',x1[-1]/3600,'2-propanol')
print('X(R)=',(rr0-rr2[-1])/rr0,'t(h)=',x2[-1]/3600,'cyclopentanol')
print('X(R)=',(rr0-rr3[-1])/rr0,'t(h)=',x3[-1]/3600,'1-oxacyclopentanol')
print('X(R)=',(rr0-rr4[-1])/rr0,'t(h)=',x4[-1]/3600,'2-oxacyclopentanol')
print('X(R)=',(rr0-rr5[-1])/rr0,'t(h)=',x5[-1]/3600,'2-adamantanol')
print('==================================')
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

