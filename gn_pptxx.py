from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdata, PESdatab
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import csv
#=====================================================================================================
# INDIVIDUAL SPECIE THERMODYNAMIC PROPERTIES COMPUTATIONS
#=====================================================================================================
#=====================================================================================================

#=====================================================================================================
# 1-OL-OXACYCLOPENTANOL
#=====================================================================================================
def ppt1(T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]

    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
    b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    adsTS_E('H2','X','HX',1,2,2,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='R' 
    print(specie_id, '  T(K)=', To)  # 1-OL-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-ol-oxacyclopentanol/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-ol-oxacyclopentanol/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='P' 
    print(specie_id, '  T(K)=', To)  # 1-OXACYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-oxacyclopentanone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-oxacyclopentanone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='RX' 
    print(specie_id, '  T(K)=', To)  # ADS FORM OF 1-OL-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/initial_adsorption/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/initial_adsorption/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='UX' 
    print(specie_id, '  T(K)=', To)  # ALKOXY ADS OF 1-OL-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    tot_mm.append(tot_mmm)
    
    specie_id='VX' 
    print(specie_id, '  T(K)=', To)  # BRIDGE FORM OF 1-OXACYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PX' 
    print(specie_id, '  T(K)=', To)  # TOP FORM OF 1-OXACYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_top/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_top/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR 1-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/TS/OH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/TS/OH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR 1-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/TS/OH-CH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/TS/OH-CH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSR'   
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF 1-OXACYCLOPENTANOL
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('R','X','RX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    """
    specie_id='TSR'   
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF 1-OXACYCLOPENTANOL
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    a = 2.701e-10 # lattice constant for Ru in \AA
    A = np.sqrt(3) * a**2 / (4 * m**2) # area of site
    adsTS_E('R','X','RX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id,tot_mm,A)
    print()
    """
    
    specie_id='TSP' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF 1-OXACYCLOPENTANONE
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('P','X','PX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print() 

    G_g1=PESdata(T,filename_link_fomat2,G_m,int_id,'X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2')  # with ads/des TS
    fig = plt.figure(1)       # <=========
    plt.plot(G_g1,label= '1-ol-oxacyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()

    fig2 = plt.figure()       # <=========
    plt.plot(G_g1,label= '1-ol-oxacyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()
    
    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   #, filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"

#=====================================================================================================
# PROPANOL
#=====================================================================================================

def ppt2(T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]
    
    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
    b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('H2','X','HX',1,2,2,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()    
    
    specie_id='Ra' 
    print(specie_id, '  T(K)=', To)  # 2-PROPANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/gas/iPrOH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/gas/iPrOH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pa' 
    print(specie_id, '  T(K)=', To)  # ACETONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/gas/acetone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/gas/acetone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='RaX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED 2-PROPANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/initial_adsorption/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/initial_adsorption/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='UaX' 
    print(specie_id, '  T(K)=', To)  # ALKOXY ADSORBATE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Alkoxy/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Alkoxy/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='VaX' 
    print(specie_id, '  T(K)=', To)  # 2-PROPANOL-BRIDGE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Ketone_bridge/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Ketone_bridge/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PaX' 
    print(specie_id, '  T(K)=', To)  # 2-PROPANOL-TOP
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Ketone_top/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/Intermediates/Ketone_top/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1a' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR 2-PROPANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/TS/OH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/TS/OH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2a' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR 2-PROPANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/TS/OH-CH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/TS/OH-CH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSRa' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF 2-PROPANOL
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('Ra','X','RaX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='TSPa' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF PROPANONE
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Pa','X','PaX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    G_g2=PESdata(T,filename_link_fomat2,G_m,int_id,'X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2')  # with ads/des TS
    fig = plt.figure(1)       # <=========
    plt.plot(G_g2,label= 'propanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()

    fig2 = plt.figure()       # <=========
    plt.plot(G_g2,label= 'propanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()
    
    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   #, filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"

#=====================================================================================================
# CYCLOPENTANOL
#=====================================================================================================
        
def ppt3(T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]
    
    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
    b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('H2','X','HX',1,2,2,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print() 
    
    specie_id='Rb' 
    print(specie_id, '  T(K)=', To)  # CYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanol/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanol/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pb' 
    print(specie_id, '  T(K)=', To)  # CYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='RbX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED CYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/initial_adsorption/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/initial_adsorption/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='UbX' 
    print(specie_id, '  T(K)=', To)  # ALKOXY ADSORBATE OF CYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Alkoxy/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Alkoxy/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='VbX' 
    print(specie_id, '  T(K)=', To)  # BRIDGED ADSORBED CYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Ketone_bridge/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Ketone_bridge/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PbX' 
    print(specie_id, '  T(K)=', To)  # TOP ADSORBED CYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Ketone_top/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/Intermediates/Ketone_top/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1b' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR THE 1ST DEH2
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/TS/OH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/TS/OH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2b' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR THE 2ND DEH2
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/TS/OH-CH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/TS/OH-CH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSRb' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF CYCLOPENTANOL
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Rb','X','RbX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='TSPb' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF CYCLOPENTANONE
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Pb','X','PbX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    G_g3=PESdata(T,filename_link_fomat2,G_m,int_id,'X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2')  # with ads/des TS
    fig = plt.figure(1)       # <=========
    plt.plot(G_g3,label= 'cyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()

    fig2 = plt.figure()       # <=========
    plt.plot(G_g3,label= 'cyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()
    
    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   #, filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"

#=====================================================================================================
# 2-OL-OXACYCLOPENTANOL
#=====================================================================================================
    
def ppt4(T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]

    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
    b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('H2','X','HX',1,2,2,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='Rc' 
    print(specie_id, '  T(K)=', To)  # 2-OXACYCLOPENTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/gas/2-ol-oxacyclopentanol/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pc' 
    print(specie_id, '  T(K)=', To)  # 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/gas/2-ol-oxacyclopentanone/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='RcX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/initial_adsorption/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='UcX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED ALKOXY OF 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='VcX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGED 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='WcX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED FROM ETHER ONLY 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_from_ether_only/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PcX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED TOP 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_top/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1c' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/TS/OH/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2c' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/TS/OH-CH/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSRc' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF 2-OXACYCLOPENTANOL
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Rc','X','RcX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='TSPc' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF 2-OXACYCLOPENTANONE
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Pc','X','PcX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    G_g4=PESdata(T,filename_link_fomat2,G_m,int_id,'X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2')  # with ads/des TS
    fig = plt.figure(1)       # <=========
    plt.plot(G_g4,label= '2-ol-oxacyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()

    fig2 = plt.figure()       # <=========
    plt.plot(G_g4,label= '2-ol-oxacyclopentanol')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()
    
    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   #, filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"

#=====================================================================================================
# 2-ADAMANTANOL
#=====================================================================================================
       
def ppt5(T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]

    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
    b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('H2','X','HX',1,2,2,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()    
    
    specie_id='Rd' 
    print(specie_id, '  T(K)=', To)  # 2-ADAMANTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/gas/adamantanol/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)  
    
    specie_id='Pd' 
    print(specie_id, '  T(K)=', To)  # 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/gas/adamantanone/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='RdX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED 2-ADAMANTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/initial_adsorption/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='UdX' 
    print(specie_id, '  T(K)=', To)  # ALKOXY FOR 2-ADAMANTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Alkoxy/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='VdX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE FCC FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/fcc/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='W1dX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE HCP-DIPOL FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/hcp/FREQ_DIPOL"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='W2dX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE HCP-NO-DIPOL FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/hcp/no-DIPOL"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PdX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED TOP FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_top/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1d' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/TS/OH/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2d' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/TS/OH-CH/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,1,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSRd' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF 2-OXACYCLOPENTANOL
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Rd','X','RdX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='TSPd' 
    print(specie_id, '  T(K)=', To)  # TSR FOR ADS OF 2-OXACYCLOPENTANONE
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('Pd','X','PdX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    G_g5=PESdata(T,filename_link_fomat2,G_m,int_id,'X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2')  # with ads/des TS
    fig = plt.figure()       # <=========
    plt.plot(G_g5,label= '2-adamantanol+Bridge+FCC')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()   

    fig2 = plt.figure()       # <=========
    plt.plot(G_g5,label= '2-adamantanol+Bridge+FCC')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show() 

    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   #, filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"



#=====================================================================================================
# SKETCHING OF 'PES-PLOT' FOR DIFFERENT CONDITIONS IN DEHYDROGENATION OF 'R' INTO 'P' & H2
#=====================================================================================================
def pptE(T):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='R' 
    print(specie_id, '  T(K)=', To)  # 1-OL-OXACYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-ol-oxacyclopentanol/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-ol-oxacyclopentanol/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='P' 
    print(specie_id, '  T(K)=', To)  # 1-OXACYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-oxacyclopentanone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/gas/1-oxacyclopentanone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
       
    #=======
    
    specie_id='Ra' 
    print(specie_id, '  T(K)=', To)  # 2-PROPANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/gas/iPrOH/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/gas/iPrOH/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pa' 
    print(specie_id, '  T(K)=', To)  # ACETONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/2-propanol/gas/acetone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/2-propanol/gas/acetone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    #========
    
    specie_id='Rb' 
    print(specie_id, '  T(K)=', To)  # CYCLOPENTANOL
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanol/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanol/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pb' 
    print(specie_id, '  T(K)=', To)  # CYCLOPENTANONE
    int_id.append(specie_id)
    a="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanone/FREQ/CONTCAR"
    b="/scratch/toyegoke/span_theory/cyclopentanol/gas/cyclopentanone/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)

    #==========
       
    specie_id='Rc' 
    print(specie_id, '  T(K)=', To)  # 2-OXACYCLOPENTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/gas/2-ol-oxacyclopentanol/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='Pc' 
    print(specie_id, '  T(K)=', To)  # 2-OXACYCLOPENTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/gas/2-ol-oxacyclopentanone/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    #===========
        
    specie_id='Rd' 
    print(specie_id, '  T(K)=', To)  # 2-ADAMANTANOL
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/gas/adamantanol/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)  
    
    specie_id='Pd' 
    print(specie_id, '  T(K)=', To)  # 2-ADAMANTANONE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/span_theory/2-adamantanol/gas/adamantanone/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
        
    return G_m, int_id, S_m, H_m
    


#=====================================================================================================
# AMMONIAL DEHYDROGENAION OVER RUTHERNIUM
#=====================================================================================================
       
def ppt6 (T, filename_link_fomat, filename_link_fomat2):
    To=T; T = To; int_id=[]; ZPE_m=[]; Elec_m=[]; S_m=[]; H_m=[]; TS_m=[]
    dH_m=[]; dH_r_m=[]; dH_t_m=[]; dH_v_m=[]; G_m=[]; s_r_m=[]; s_t_m=[]
    s_v_m=[]; ts_r_m=[]; ts_t_m=[]; ts_v_m=[]; dH_TS_m=[]; tot_mm=[]; q_t_m=[]
    q_r_m=[]; q_v_m=[]

    specie_id='X' 
    print(specie_id, '  T(K)=', To)  # Ru CATALYST
    int_id.append(specie_id)
    a="/scratch/toyegoke/NH3_to_N_pathway/bare-slab/CONTCAR"
    b="/scratch/toyegoke/NH3_to_N_pathway/bare-slab/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='HX' 
    print(specie_id, '  T(K)=', To)  # ADSROBED HYDROGEN ATOM
    int_id.append(specie_id)
    a="/scratch/toyegoke/NH3_to_N_pathway/1H-Ru0001/FREQ/CONTCAR"
    b="/scratch/toyegoke/NH3_to_N_pathway/1H-Ru0001/FREQ/OUTCAR"
    ccc="/scratch/toyegoke/NH3_to_N_pathway/1H-Ru0001/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id,ccc)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='H2' 
    print(specie_id, '  T(K)=', To)  # HYDROGEN GAS
    int_id.append(specie_id)
    a="/scratch/toyegoke/NH3_to_N_pathway/H2-gas/FREQ/CONTCAR"
    b="/scratch/toyegoke/NH3_to_N_pathway/H2-gas/FREQ/OUTCAR"
    # a="/scratch/toyegoke/span_theory/H2/FREQ/CONTCAR"
    # b="/scratch/toyegoke/span_theory/H2/FREQ/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSH2' 
    print(specie_id, '  T(K)=', To)   # TSR FOR ADS OF HYDROGEN GAS
    int_id.append(specie_id)
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    adsTS_E('H2','X','HX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()    
    
    specie_id='R' 
    print(specie_id, '  T(K)=', To)  # AMMONIA
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH3-gas/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)  
    
    specie_id='P' 
    print(specie_id, '  T(K)=', To)  # NITROGEN-GAS
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/N2-gas/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='RX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED AMMONIA
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH3/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='UX' 
    print(specie_id, '  T(K)=', To)  # NH2-INTERMEDIATE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH2/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='VX' 
    print(specie_id, '  T(K)=', To)  # NH-INTERMEDIATE
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='PX' 
    print(specie_id, '  T(K)=', To)  # ADSORBED NITROGEN
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/N/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS1' 
    print(specie_id, '  T(K)=', To)  # TS1 FOR NH3
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH2-H_TS/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS2' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR NH2
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/NH-H_TS/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TS3' 
    print(specie_id, '  T(K)=', To)  # TS2 FOR NH
    int_id.append(specie_id)
    xy="/scratch/toyegoke/NH3_to_N_pathway/N-H_TS/FREQ"
    a=xy+"/CONTCAR"
    b=xy+"/OUTCAR"
    [S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,3,To,specie_id)
    S_m.append(S)
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    H_m.append(H)
    G_m.append(G)
    dH_t_m.append(dH_t)
    
    specie_id='TSR' 
    print(specie_id, '  T(K)=', To)  # TSR FOR NH3
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('R','X','RX',2,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    specie_id='TSP' 
    print(specie_id, '  T(K)=', To)  # TSR FOR N2
    s_t=0
    dH_t=0
    s_v=0
    dH_v=0
    s_v_m.append(s_v)
    dH_v_m.append(dH_v)
    s_t_m.append(s_t)
    dH_t_m.append(dH_t)
    int_id.append(specie_id)
    adsTS_E('P','X','PX',1,1,1,T,G_m, H_m, dH_t_m, S_m, s_t_m,int_id)
    print()
    
    G_g6=PESdatab(T,filename_link_fomat2,G_m,int_id,'X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TS3','TSR','TSP','TSH2')  # with ads/des TS
    fig = plt.figure()       # <=========
    plt.plot(G_g6,label= 'Ammonia Dehydrogenation over Ru')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show()   

    fig2 = plt.figure()       # <=========
    plt.plot(G_g6,label= 'Ammonia Dehydrogenation over Ru')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlabel('Reaction Path')
    plt.ylabel('Energy in eV')
    plt.xticks(rotation=90)
    plt.title('Energy Profile at T(K) = {}'.format(round(T,2)))
    plt.legend()
    #plt.show() 

    csvfile = filename_link_fomat #"/home/toyegoke/ENS_kinetics/"+title+"_THERMO_data_.csv"
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),'G_x(J/mol) ','H_x (J/mol)','S_x (J/mol/K)','G_x(eV) ','H_x (eV)','S_x (eV/K)'])
        glen=len(G_m)
        for i in range(glen):            
            wr.writerow([int_id[i],G_m[i],H_m[i],S_m[i],G_m[i]/96485,H_m[i]/96485,S_m[i]/96485])

    return G_m, int_id, S_m, H_m, fig, fig2   

