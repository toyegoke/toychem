from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES 
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const

#=====================================================================================================
#=====================================================================================================
# CONDITIONS FOR THE THERMODYNMAIC PPTY COMPUTATIONS // INPUT DETAILS // DATA COLLECTIONS
#=====================================================================================================
#=====================================================================================================


To=550 #273.15+25; 
T = To

int_id=[]
ZPE_m=[]
Elec_m=[]
S_m=[]
H_m=[]
TS_m=[]
dH_m=[]
dH_r_m=[]
dH_t_m=[]
dH_v_m=[]
G_m=[]
s_r_m=[]
s_t_m=[]
s_v_m=[]
ts_r_m=[]
ts_t_m=[]
ts_v_m=[]
dH_TS_m=[]
tot_mm=[]
q_t_m=[]
q_r_m=[] 
q_v_m=[]


#=====================================================================================================
#=====================================================================================================
# INDIVIDUAL SPECIE THERMODYNAMIC PROPERTIES COMPUTATIONS
#=====================================================================================================
#=====================================================================================================


specie_id='X' 
print(specie_id, '  T(K)=', To)  # Ru CATALYST
int_id.append(specie_id)
a="/home/toyegoke/span_theory/bare-slab/CONTCAR"
b="/home/toyegoke/span_theory/bare-slab/OUTCAR"
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
a="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
ccc="/scratch/toyegoke/span_theory/1H_Ru0001/FREQ/OUTCAR"
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


#=====================================================================================================
#=====================================================================================================

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


specie_id='RaX' 
print(specie_id, '  T(K)=', To)  # ADSORBED 2-PROPANOL
int_id.append(specie_id)
a="/scratch/toyegoke/span_theory/2-propanol/initial_adsorption/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/2-propanol/initial_adsorption/FREQ/OUTCAR"
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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



#=====================================================================================================
#=====================================================================================================

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


specie_id='RdX' 
print(specie_id, '  T(K)=', To)  # ADSORBED 2-ADAMANTANOL
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/initial_adsorption/FREQ"
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


specie_id='UdX' 
print(specie_id, '  T(K)=', To)  # ALKOXY FOR 2-ADAMANTANOL
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Alkoxy/FREQ"
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


specie_id='VdX' 
print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE FCC FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/fcc/FREQ"
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


specie_id='W1dX' 
print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE HCP-DIPOL FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/hcp/FREQ_DIPOL"
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


specie_id='W2dX' 
print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGE HCP-NO-DIPOL FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_bridge/hcp/no-DIPOL"
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


specie_id='PdX' 
print(specie_id, '  T(K)=', To)  # ADSORBED TOP FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/Intermediates/Ketone_top/FREQ"
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


specie_id='TS1d' 
print(specie_id, '  T(K)=', To)  # TS1 FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/TS/OH/FREQ"
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


specie_id='TS2d' 
print(specie_id, '  T(K)=', To)  # TS2 FOR 2-ADAMANTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-adamantanol/TS/OH-CH/FREQ"
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


#=====================================================================================================
#=====================================================================================================

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


specie_id='RcX' 
print(specie_id, '  T(K)=', To)  # ADSORBED 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/initial_adsorption/FREQ"
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


specie_id='UcX' 
print(specie_id, '  T(K)=', To)  # ADSORBED ALKOXY OF 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ"
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


specie_id='VcX' 
print(specie_id, '  T(K)=', To)  # ADSORBED BRIDGED 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ"
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


specie_id='WcX' 
print(specie_id, '  T(K)=', To)  # ADSORBED FROM ETHER ONLY 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_from_ether_only/FREQ"
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


specie_id='PcX' 
print(specie_id, '  T(K)=', To)  # ADSORBED TOP 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/Intermediates/Ketone_top/FREQ"
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


specie_id='TS1c' 
print(specie_id, '  T(K)=', To)  # TS1 FOR 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/TS/OH/FREQ"
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


specie_id='TS2c' 
print(specie_id, '  T(K)=', To)  # TS2 FOR 2-OXACYCLOPENTANONE
int_id.append(specie_id)
xy="/scratch/toyegoke/span_theory/2-ol-oxacyclopentanol/TS/OH-CH/FREQ"
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

#=====================================================================================================
#=====================================================================================================

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


specie_id='RbX' 
print(specie_id, '  T(K)=', To)  # ADSORBED CYCLOPENTANOL
int_id.append(specie_id)
a="/scratch/toyegoke/span_theory/cyclopentanol/initial_adsorption/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/cyclopentanol/initial_adsorption/FREQ/OUTCAR"
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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




#=====================================================================================================
#=====================================================================================================


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


specie_id='RX' 
print(specie_id, '  T(K)=', To)  # ADS FORM OF 1-OL-OXACYCLOPENTANOL
int_id.append(specie_id)
a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/initial_adsorption/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/initial_adsorption/FREQ/OUTCAR"
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
S_m.append(S)
s_v_m.append(s_v)
dH_v_m.append(dH_v)
s_t_m.append(s_t)
H_m.append(H)
G_m.append(G)
dH_t_m.append(dH_t)


specie_id='UX' 
print(specie_id, '  T(K)=', To)  # ALKOXY ADS OF 1-OL-OXACYCLOPENTANOL
int_id.append(specie_id)
a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Alkoxy/FREQ/OUTCAR"
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
S_m.append(S)
s_v_m.append(s_v)
dH_v_m.append(dH_v)
s_t_m.append(s_t)
H_m.append(H)
G_m.append(G)
dH_t_m.append(dH_t)


specie_id='VX' 
print(specie_id, '  T(K)=', To)  # BRIDGE FORM OF 1-OXACYCLOPENTANONE
int_id.append(specie_id)
a="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ/CONTCAR"
b="/scratch/toyegoke/span_theory/1-ol-oxacyclopentanol/Intermediates/Ketone_bridge/FREQ/OUTCAR"
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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
[S, Cp, dH,  Q,  G, tetha_rot, tot_mmm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v]=thermoppt(a,b,2,To,specie_id)
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




#=====================================================================================================
#=====================================================================================================
# SKETCHING OF 'PES-PLOT' FOR DIFFERENT CONDITIONS IN DEHYDROGENATION OF 'R' INTO 'P' & H2
#=====================================================================================================
#=====================================================================================================

# FREE ENERGY (G)
G_g1=PES(G_m,int_id,'X','R','RX','UX','HX','VX','PX','H2','P','TS1','TS2','TSR','TSP','TSH2')  # with ads/des TS
G_g2=PES(G_m,int_id,'X','Ra','RaX','UaX','HX','VaX','PaX','H2','Pa','TS1a','TS2a','TSRa','TSPa','TSH2')  # with ads/des TS
G_g3=PES(G_m,int_id,'X','Rb','RbX','UbX','HX','VbX','PbX','H2','Pb','TS1b','TS2b','TSRb','TSPb','TSH2')  # with ads/des TS
G_g4=PES(G_m,int_id,'X','Rc','RcX','UcX','HX','VcX','PcX','H2','Pc','TS1c','TS2c','TSRc','TSPc','TSH2')  # with ads/des TS
G_g5=PES(G_m,int_id,'X','Rd','RdX','UdX','HX','VdX','PdX','H2','Pd','TS1d','TS2d','TSRd','TSPd','TSH2')  # with ads/des TS
   
plt.figure(1)
plt.plot(G_g2,label= 'propanol')
plt.plot(G_g3,label= 'cyclopentanol')
plt.plot(G_g1,label= '1-ol-oxacyclopentanol')
plt.plot(G_g4,label= '2-ol-oxacyclopentanol')
plt.plot(G_g5,label= '2-adamantanol+Bridge+FCC')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.xlabel('Reaction Path')
plt.ylabel('Energy in eV')
plt.xticks(rotation=90)
plt.legend()
#plt.show()