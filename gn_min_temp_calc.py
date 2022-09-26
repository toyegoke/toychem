from gn_thermop import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES 
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5, pptE
####from kin_temp0 import kineticparameter, dydt,solve_odes
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import csv
from numpy import round

from gn_thermop import specie_index

# GOOD FOR ESTIMATING MINIMUM THERMO CONDITION (dG, dH, dS, TEMP) FOR REACTION IS POSSIBLE FOR DIFF ALCOHOLS
#ergornicity
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


def temp_effect(Order_specieID_list, G_m0):      # <=============  

    TT = np.linspace(0.0001, 600, 10) 

    ppt = pptE
    Eox1 = []
    Epro = []
    Ecyc = []
    Eox2 = []
    Eadm = []
    S_x_in = (18*1e-6)# mol/m^2 (equivalent of 0.44 ML from expt data)
    S_x = S_x_in/(1e-6) # umol/m^2 (equivalent of 0.44 ML from expt data)
    tlens=len(TT)
    i=0
    if G_m0=='G_m':
        XX_0 = 0
    elif G_m0=='S_m':
        XX_0 = 1
    elif G_m0=='H_m':
        XX_0 = 2

    for i in range(tlens):

        T=TT[i]
        G_m, int_id, S_m, H_m = ppt(T)      # <============= 
        XX_1 = [G_m, S_m, H_m]
        energylist=XX_1[XX_0]   
                
        # list of selected specie IDs in a special order with string  elements
        rxn_specie_id = Order_specieID_list
        h2=rxn_specie_id[0]
        rr=rxn_specie_id[1] 
        pp=rxn_specie_id[2] 
        ra=rxn_specie_id[3]
        pa=rxn_specie_id[4]
        rb=rxn_specie_id[5]
        pb=rxn_specie_id[6] 
        rc=rxn_specie_id[7]
        pc=rxn_specie_id[8]
        rd=rxn_specie_id[9]
        pd=rxn_specie_id[10]
        print(h2,rr,pp,ra,pa,rb,pb,rc,pc,rd,pd)
        print(int_id)
        specie_index(h2,int_id)
        specie_index(rr,int_id)
        specie_index(pp,int_id)
        specie_index(ra,int_id)
        specie_index(pa,int_id)
        print(G_m[specie_index(ra,int_id)])
        
        #Free Gibb (G) energy for the species
        H2g=energylist[specie_index(h2,int_id)]
        Rgx =energylist[specie_index(rr,int_id)]
        Pgx =energylist[specie_index(pp,int_id)]
        Rga=energylist[specie_index('Ra',int_id)]
        Pga=energylist[specie_index(pa,int_id)]
        Rgb=energylist[specie_index(rb,int_id)]
        Pgb=energylist[specie_index(pb,int_id)]
        Rgc=energylist[specie_index(rc,int_id)]
        Pgc=energylist[specie_index(pc,int_id)]
        Rgd=energylist[specie_index(rd,int_id)]
        Pgd=energylist[specie_index(pd,int_id)]
        
        #Free Gibb (G) energy for the species
        RxnXx=Pgx+H2g-Rgx
        RxnXa=Pga+H2g-Rga
        RxnXb=Pgb+H2g-Rgb
        RxnXc=Pgc+H2g-Rgc
        RxnXd=Pgd+H2g-Rgd
        Eox1.append(RxnXx/96485)
        Epro.append(RxnXa/96485)
        Ecyc.append(RxnXb/96485)
        Eox2.append(RxnXc/96485)
        Eadm.append(RxnXd/96485)
	
    return TT,Eox1,Epro,Ecyc,Eox2,Eadm    


title2= 'propanol'
title3= 'cyclopentanol'
title1= '1-ol-oxacyclopentanol'
title4= '2-ol-oxacyclopentanol'
title5= '2-adamantanol+Bridge+FCC' 

Exx = 'G_m'
TT,Eox1,Epro,Ecyc,Eox2,Eadm=temp_effect(['H2','R','P','Ra','Pa','Rb','Pb','Rc','Pc','Rd','Pd'], Exx)
slope1, intercept1 = np.polyfit(TT, Eox1, 1)
slope2, intercept2 = np.polyfit(TT, Epro, 1)
slope3, intercept3 = np.polyfit(TT, Ecyc, 1)
slope4, intercept4 = np.polyfit(TT, Eox2, 1)
slope5, intercept5 = np.polyfit(TT, Eadm, 1)
slope=['slope',slope1,slope2,slope3,slope4,slope5]
intercept=['intercept',intercept1,intercept2,intercept3,intercept4,intercept5]
ppp = [TT,Eox1,Epro,Ecyc,Eox2,Eadm]
ranP = len(TT)
Aaa=Exx.split('_')
qua=Aaa[0]
csvfile = "/home/toyegoke/ENS_kinetics/Final-min-temp/data_"+qua+"_for_species.csv"
with open(csvfile, "a") as fp:
    wr = csv.writer(fp, dialect='excel')
    wr.writerow(['Temp','Eox1_'+qua,'Epro_'+qua,'Ecyc_'+qua,'Eox2_'+qua,'Eadm_'+qua])
    wr.writerow([TT[0],Eox1[0],Epro[0],Ecyc[0],Eox2[0],Eadm[0]])
    wr.writerow([TT[1],Eox1[1],Epro[1],Ecyc[1],Eox2[1],Eadm[1]])
    wr.writerow([TT[2],Eox1[2],Epro[2],Ecyc[2],Eox2[2],Eadm[2]])
    wr.writerow([TT[3],Eox1[3],Epro[3],Ecyc[3],Eox2[3],Eadm[3]])
    wr.writerow([TT[4],Eox1[4],Epro[4],Ecyc[4],Eox2[4],Eadm[4]])
    wr.writerow([TT[5],Eox1[5],Epro[5],Ecyc[5],Eox2[5],Eadm[5]])
    wr.writerow([TT[6],Eox1[6],Epro[6],Ecyc[6],Eox2[6],Eadm[6]])
    wr.writerow([TT[7],Eox1[7],Epro[7],Ecyc[7],Eox2[7],Eadm[7]])
    wr.writerow([TT[8],Eox1[8],Epro[8],Ecyc[8],Eox2[8],Eadm[8]])
    wr.writerow([TT[9],Eox1[9],Epro[9],Ecyc[9],Eox2[9],Eadm[9]])
    wr.writerow([slope[0],slope[1],slope[2],slope[3],slope[4],slope[5]])
    wr.writerow([intercept[0],intercept[1],intercept[2],intercept[3],intercept[4],intercept[5]])
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(TT, Eox1, linestyle='--', label='1-ol-oxacyclopentanol')
ax.plot(TT, Epro, linestyle='--', label='Propanol')
ax.plot(TT, Ecyc, linestyle='--', label='Cyclopentanol')
ax.plot(TT, Eox2, linestyle='--', label='2-ol-oxacyclopentanol')
ax.plot(TT, Eadm, linestyle='--', label='2-adamantanol')
ax.plot(TT, np.zeros(ranP), 'k')
plt.xlabel('T / K')
plt.ylabel('$\Delta$'+qua+' / eV')
#plt.title('Free Energy')  
ax.legend()
plt.show()
fig.savefig('/home/toyegoke/ENS_kinetics/Final-min-temp/plotdata_'+qua+'_for_species.png')

Exx = 'H_m'
TT,Eox1,Epro,Ecyc,Eox2,Eadm=temp_effect(['H2','R','P','Ra','Pa','Rb','Pb','Rc','Pc','Rd','Pd'], Exx)
slope1, intercept1 = np.polyfit(TT, Eox1, 1)
slope2, intercept2 = np.polyfit(TT, Epro, 1)
slope3, intercept3 = np.polyfit(TT, Ecyc, 1)
slope4, intercept4 = np.polyfit(TT, Eox2, 1)
slope5, intercept5 = np.polyfit(TT, Eadm, 1)
slope=['slope',slope1,slope2,slope3,slope4,slope5]
intercept=['intercept',intercept1,intercept2,intercept3,intercept4,intercept5]
ppp = [TT,Eox1,Epro,Ecyc,Eox2,Eadm]
ranP = len(TT)
Aaa=Exx.split('_')
qua=Aaa[0]
csvfile = "/home/toyegoke/ENS_kinetics/Final-min-temp/data_"+qua+"_for_species.csv"
with open(csvfile, "a") as fp:
    wr = csv.writer(fp, dialect='excel')
    wr.writerow(['Temp','Eox1_'+qua,'Epro_'+qua,'Ecyc_'+qua,'Eox2_'+qua,'Eadm_'+qua])
    wr.writerow([TT[0],Eox1[0],Epro[0],Ecyc[0],Eox2[0],Eadm[0]])
    wr.writerow([TT[1],Eox1[1],Epro[1],Ecyc[1],Eox2[1],Eadm[1]])
    wr.writerow([TT[2],Eox1[2],Epro[2],Ecyc[2],Eox2[2],Eadm[2]])
    wr.writerow([TT[3],Eox1[3],Epro[3],Ecyc[3],Eox2[3],Eadm[3]])
    wr.writerow([TT[4],Eox1[4],Epro[4],Ecyc[4],Eox2[4],Eadm[4]])
    wr.writerow([TT[5],Eox1[5],Epro[5],Ecyc[5],Eox2[5],Eadm[5]])
    wr.writerow([TT[6],Eox1[6],Epro[6],Ecyc[6],Eox2[6],Eadm[6]])
    wr.writerow([TT[7],Eox1[7],Epro[7],Ecyc[7],Eox2[7],Eadm[7]])
    wr.writerow([TT[8],Eox1[8],Epro[8],Ecyc[8],Eox2[8],Eadm[8]])
    wr.writerow([TT[9],Eox1[9],Epro[9],Ecyc[9],Eox2[9],Eadm[9]])
    wr.writerow([slope[0],slope[1],slope[2],slope[3],slope[4],slope[5]])
    wr.writerow([intercept[0],intercept[1],intercept[2],intercept[3],intercept[4],intercept[5]])
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(TT, Eox1, linestyle='--', label='1-ol-oxacyclopentanol')
ax.plot(TT, Epro, linestyle='--', label='Propanol')
ax.plot(TT, Ecyc, linestyle='--', label='Cyclopentanol')
ax.plot(TT, Eox2, linestyle='--', label='2-ol-oxacyclopentanol')
ax.plot(TT, Eadm, linestyle='--', label='2-adamantanol')
ax.plot(TT, np.zeros(ranP), 'k')
plt.xlabel('T / K')
plt.ylabel('$\Delta$'+qua+' / eV')
#plt.title('Enthalpy')
ax.legend()     
plt.show()
fig.savefig('/home/toyegoke/ENS_kinetics/Final-min-temp/plotdata_'+qua+'_for_species.png')


Exx = 'S_m'
TT,Eox1,Epro,Ecyc,Eox2,Eadm=temp_effect(['H2','R','P','Ra','Pa','Rb','Pb','Rc','Pc','Rd','Pd'], Exx)
slope1, intercept1 = np.polyfit(TT, Eox1, 1)
slope2, intercept2 = np.polyfit(TT, Epro, 1)
slope3, intercept3 = np.polyfit(TT, Ecyc, 1)
slope4, intercept4 = np.polyfit(TT, Eox2, 1)
slope5, intercept5 = np.polyfit(TT, Eadm, 1)
slope=['slope',slope1,slope2,slope3,slope4,slope5]
intercept=['intercept',intercept1,intercept2,intercept3,intercept4,intercept5]
ppp = [TT,Eox1,Epro,Ecyc,Eox2,Eadm]
ranP = len(TT)
Aaa=Exx.split('_')
qua=Aaa[0]
csvfile = "/home/toyegoke/ENS_kinetics/Final-min-temp/data_"+qua+"_for_species.csv"
with open(csvfile, "a") as fp:
    wr = csv.writer(fp, dialect='excel')
    wr.writerow(['Temp','Eox1_'+qua,'Epro_'+qua,'Ecyc_'+qua,'Eox2_'+qua,'Eadm_'+qua])
    wr.writerow([TT[0],Eox1[0],Epro[0],Ecyc[0],Eox2[0],Eadm[0]])
    wr.writerow([TT[1],Eox1[1],Epro[1],Ecyc[1],Eox2[1],Eadm[1]])
    wr.writerow([TT[2],Eox1[2],Epro[2],Ecyc[2],Eox2[2],Eadm[2]])
    wr.writerow([TT[3],Eox1[3],Epro[3],Ecyc[3],Eox2[3],Eadm[3]])
    wr.writerow([TT[4],Eox1[4],Epro[4],Ecyc[4],Eox2[4],Eadm[4]])
    wr.writerow([TT[5],Eox1[5],Epro[5],Ecyc[5],Eox2[5],Eadm[5]])
    wr.writerow([TT[6],Eox1[6],Epro[6],Ecyc[6],Eox2[6],Eadm[6]])
    wr.writerow([TT[7],Eox1[7],Epro[7],Ecyc[7],Eox2[7],Eadm[7]])
    wr.writerow([TT[8],Eox1[8],Epro[8],Ecyc[8],Eox2[8],Eadm[8]])
    wr.writerow([TT[9],Eox1[9],Epro[9],Ecyc[9],Eox2[9],Eadm[9]])
    wr.writerow([slope[0],slope[1],slope[2],slope[3],slope[4],slope[5]])
    wr.writerow([intercept[0],intercept[1],intercept[2],intercept[3],intercept[4],intercept[5]])
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(TT, Eox1, linestyle='--', label='1-ol-oxacyclopentanol')
ax.plot(TT, Epro, linestyle='--', label='Propanol')
ax.plot(TT, Ecyc, linestyle='--', label='Cyclopentanol')
ax.plot(TT, Eox2, linestyle='--', label='2-ol-oxacyclopentanol')
ax.plot(TT, Eadm, linestyle='--', label='2-adamantanol')
ax.plot(TT, np.zeros(ranP), 'k')
plt.xlabel('T / K')
plt.ylabel('$\Delta$'+qua+' / eV')
#plt.title('Entropy')
ax.legend()
plt.show()
fig.savefig('/home/toyegoke/ENS_kinetics/Final-min-temp/plotdata_'+qua+'_for_species.png')

print('COMPLETED')




