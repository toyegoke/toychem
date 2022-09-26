# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:22:05 2021

@author: Toyese OYEGOKE, Stephane STEINMANN, Carine MICHEL
@affiliation: Laboratoire de Chimie, ENS Lyon, FRANCE
"""

from gn_thermop import thermoppt 
#from gn_thermop_copy import thermoppt 
from gn_thermop import adsTS_E 
from gn_thermop import PES, PESdata, ESPAN_PES
from gn_pptx import ppt1,ppt2,ppt3,ppt4,ppt5
## from kin_temp0 import kineticparameter, dydt,solve_odes,solve_odes2,solve_odes3,solve_odes4,solve_odes46 #40
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
from numpy import round
from statistics import stdev, mean

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# ENERGY SPAN CALCULATION
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def espan_calc(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout):
    
    #for specific time range GOOD FOR RDS (for RDS when rate is measured in terms of Rate_x)
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime # (14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp #550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = cmass #40  # in g
    C_R0 = 997.74 # mol/m3 
    outflow_rate = flowout * (14.27) * 1e-7 # in m3/sec
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
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)
    
    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)
    
    fig.savefig('/home/toyegoke/ENS_kinetics/Final-espan/'+title+'_'+str(temp)+'K'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-espan/'+title+'_'+str(temp)+'K'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T) 

    rxn_specie_id = Order_specieID_list     # list of selected specie IDs in a special order with string  elements
    x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
    vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
    tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13]       
    
    Gs = PESdata(temp,filename_link_fomat,G_m,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2) # in eV (for the output)
    G_g = (Gs[0],Gs[2],Gs[4],Gs[6],Gs[8],Gs[10],
           Gs[12],max(Gs[12],Gs[14]),Gs[14],Gs[16],Gs[18],Gs[20],Gs[22])
    StepID = ['R+3X','TSR+2X','RX+2X','TS1+2X','UX+HX+X','TS2+HX+X','VX+2HX',
              'TSV+2HX','PX+2HX','TSP+2HX','P+2HX+X','TSH2+P+X','P+H2+3X']
    
    print('G_g=',G_g)
    print('Gs=',Gs)
    prod = G_g[-1]
    rxt = G_g[0]
    Grxn = prod - rxt
    Glen = len(G_g)
    i = 0
    Espanlist=[]
    EspanISlist=[]
    EspanTSlist=[]
    
    TSid = [1,3,5,7,9,11]
    ISid = [0,2,4,6,8,10,12]
    tslen = len(TSid)
    islen = len(ISid)
    
    j=0
    i=0
    
    for i in range(tslen):
        ts = TSid[i]# 1+2*i
        TS = G_g[ts]
        for j in range(islen):
            imtd = ISid[j]
            if imtd < ts:
                print('before TS point')
                if imtd < (Glen):                 
          
                    Gis = G_g[imtd]
                    Espan0 = TS-Gis
                    Espan = Espan0
                    Espanlist.append(Espan)
                    EspanISlist.append(imtd)
                    EspanTSlist.append(ts)
                    print('==================================')
                    print('TS=',TS)
                    print('Gis=',Gis)
                    print('ts=',ts)
                    print('imtd=',imtd)
                    print('Espan=', Espan)
                else:
                    print('outside the provided range')            
            else:#if imtd > ts:
                print('after TS point')
                if imtd < (Glen): 
                    Gis = G_g[imtd]
                    Espan0 = TS-Gis+Grxn
                    Espan = Espan0
                    Espanlist.append(Espan)
                    EspanISlist.append(imtd)
                    EspanTSlist.append(ts)
                    print('==================================')
                    print('TS=',TS)
                    print('Gis=',Gis)
                    print('ts=',ts)
                    print('imtd=',imtd)
                    print('Espan=', Espan)
                else:
                    print('outside the provided range')
    print('Espanlist=',Espanlist)
    print('EspanTSlist=',EspanTSlist)
    print('EspanISlist=',EspanISlist)
    Emax=max(Espanlist)
    Emax_id=Espanlist.index(Emax)
    Espan = Emax
    EspanISid=EspanISlist[Emax_id]                  
    EspanTSid=EspanTSlist[Emax_id]
    EspanISv=G_g[EspanISid]                  
    EspanTSv=G_g[EspanTSid]
    print('Emax=',Emax)
    print('Emax_id=',Emax_id)
    print('EspanTSid=',EspanTSid)
    print('EspanISid=',EspanISid)
    print('Espan/1000, EspanISid, EspanTSid, StepID=',Espan/1000, EspanISid, EspanTSid, StepID)
    return Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn   # to make it kJ/mol

"""
#================================================================================================
# SOLUTIONS 1
#================================================================================================
kB=const.Boltzmann # J/K
N=const.Avogadro # mol^-1
h=const.Planck # J.sam nie-7
R=N*kB # const.gas_constant in J/mol/K    

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

Estep = ['R','TSR','RX','TS2','UX+HX']

pt = [ppt1,ppt2,ppt3,ppt4,ppt5]
titl = [title1,title2,title3,title4,title5]
specie_ID = [engg1,engg2,engg3,engg4,engg5]

cmass = 1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

rtemp = [530,540,550,560,570] # 550K
rtime = 10**(5) #for 5%/20%

ranP = len(titl)
ranT = len(rtemp)
i=0


for i in range(ranP):
    ppt=pt[i]
    title=titl[i]
    engg=specie_ID[i]
    j=0
    AA=[]; BB=[]; CC=[]; DD=[]; EE=[]; FF=[]; GG=[]; TS=[]; IS=[]; RR=[]
     
    for j in range(ranT):
        
        TEMP = rtemp[j]
    
        Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn = espan_calc(ppt,title,engg,rtime,TEMP,cmass,flowout)
        
        AA.append(Espan*96.485)
        BB.append(EspanISid)
        CC.append(EspanTSid)
        DD.append(TEMP)
        EE.append(StepID[EspanISid])
        FF.append(StepID[EspanTSid])
        GG.append(Grxn*96.485)
        IS.append(EspanISv*96.485)
        TS.append(EspanTSv*96.485)
        RR.append(kB*TEMP/h*np.exp(-Espan*96485/R/TEMP))

    Gstd = stdev(GG)

    Astd = stdev(AA)
    Amx = mean(AA)    

    Rstd = stdev(RR)
    Rmx = mean(RR)   

    csvfile = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+"_data_espan_.csv"  
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['E-span(kJ/mol)', 'TS_step', 'TDTS(kJ/mol)', 'TD-TS', 'IS_step', 'TDIS(kJ/mol)', 'TD-IS', 'Grxn(kJ/mol)', 'Temp in K', 'Rate(TOF)'])
        #glen=len(ZP_g)
        for j in range(ranT):            
            wr.writerow([ AA[j], CC[j], TS[j], FF[j], BB[j], IS[j], EE[j], GG[j], DD[j], RR[j] ])
        wr.writerow(['E-span-Stdev', Astd, 'E-span-Mean', Amx, 'Grxn', Gstd, 'Rate(TOF)', Rstd])

    slope, intercept = np.polyfit(DD, AA, 1)
    mx = round(slope, 2, out=None)
    cy = round(intercept, 2, out=None)
        
    fig = plt.figure()
    ay = plt.subplot(111)
    ay.semilogx(DD, AA, label=title+' (m='+str(mx)+')')
    plt.title('Energy Span Profile')
    plt.xlabel('Temperature, $ T/K $')
    plt.ylabel('E-span, $ kJ/mol $')   
    ay.legend()  
    fig.savefig('/home/toyegoke/ENS_kinetics/Final-espan/'+title+'_E-span_plot_.png')
#================================================================================================
"""
  
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# TOF CALCULATION
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def xtof_calc(moleculer_property_function, title_of_molecule, Order_specieID_list, rtime, rtemp, cmass, flowout):
    
    #for specific time range GOOD FOR RDS (for RDS when rate is measured in terms of Rate_x)
    
    title0 = title_of_molecule  # name of the molecule studied
    TIME = rtime # (14.49) *60 # in sec # selected_time # MAXtime in hour, no of steps when using zero as low take 0.001 as 0
    temp = rtemp #550 # we choose that but optim result is 523.34 # in K # temperature_selected # for each molecule study
    mc0 = cmass #40  # in g
    C_R0 = 997.74 # mol/m3 
    outflow_rate = flowout * (14.27) * 1e-7 # in m3/sec
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
    
    filename_link_fomat = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+'_'+"_THERMO_data_.csv"     #+str(T)
    filename_link_fomat2 = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+'_'+"_PES_profile_data_.csv"  #+str(T)
    
    G_m, int_id, S_m, H_m, fig, fig2 = ppt(temp, filename_link_fomat, filename_link_fomat2)
    
    fig.savefig('/home/toyegoke/ENS_kinetics/Final-espan/'+title+'_'+str(temp)+'K'+'_plot_data_for_COMBINED_ENERGY_PROFILE.png')  # +str(T)
    fig2.savefig('/home/toyegoke/ENS_kinetics/Final-espan/'+title+'_'+str(temp)+'K'+'_plot_data_for_ENERGY_PROFILE.png')  #  +str(T) 

    rxn_specie_id = Order_specieID_list     # list of selected specie IDs in a special order with string  elements
    x=rxn_specie_id[0]; r=rxn_specie_id[1]; rx=rxn_specie_id[2]; ux=rxn_specie_id[3]; hx=rxn_specie_id[4] 
    vx=rxn_specie_id[5]; px=rxn_specie_id[6]; h2=rxn_specie_id[7]; p=rxn_specie_id[8]; ts1=rxn_specie_id[9]; ts2=rxn_specie_id[10]
    tsr=rxn_specie_id[11]; tsp=rxn_specie_id[12]; tsh2=rxn_specie_id[13]       
    
    Gs = PESdata(temp,filename_link_fomat,G_m,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2)  # in eV (from the solution)
    G_g = (Gs[0],Gs[2],Gs[4],Gs[6],Gs[8],Gs[10],
           Gs[12],max(Gs[12],Gs[14]),Gs[14],Gs[16],Gs[18],Gs[20],Gs[22])
    StepID = ['R+3X','TSR+2X','RX+2X','TS1+2X','UX+HX+X','TS2+HX+X','VX+2HX',
          'TSV+2HX','PX+2HX','TSP+2HX','P+2HX+X','TSH2+P+X','P+H2+3X']
    #StepTS = ['TSR+2X','TS1+2X','TS2+HX+X','TSV+2HX','TSP+2HX','TSH2+P+X']
    #StepIS = ['R+3X','RX+2X','UX+HX+X','VX+2HX','PX+2HX','P+2HX+X','P+H2+3X']
    StepTS = []
    StepIS = []
    
    Stepis = []
    Stepts = []
    Steptemp = []
    Steptof = []
    StepEis = []
    StepEts = []
    StepE = []
    
    print('G_g=',G_g);     print('Gs=',Gs)
    prod = G_g[-1];     rxt = G_g[0]
    Grxn = prod - rxt
    Glen = len(G_g)
    Espanlist=[]
    EspanISlist=[]
    EspanTSlist=[]
    
    TSid = [1,3,5,7,9,11]

    ISid = [0,2,4,6,8,10,12] # when product and reactant step are included
    #ISid = [2,4,6,8,10] # when product and reactant step are NOT included    
    
    tslen = len(TSid)
    islen = len(ISid)  

    j=0
    i=0
    TSTOF=[]
    
    for i in range(tslen):
        ts = TSid[i]# 1+2*i
        TS = G_g[ts]
        RR = []; rr = []
        for j in range(islen):
            imtd = ISid[j]
            if imtd > ts:
                print('after TS point')
                if imtd < (Glen):                 
                    Gis = G_g[imtd]
                    Espan0 = TS-Gis-abs(Grxn)
                    Espan0 = TS-Gis+Grxn
                    #Espan1 = TS-Gis   # the energy in eV
                    Espan = Espan0
                    Espanlist.append(Espan)
                    EspanISlist.append(imtd)
                    EspanTSlist.append(ts)
                    StepIS.append(StepID[imtd])
                    StepTS.append(StepID[ts])   
                    RR.append(np.exp(Espan1*96485/R/temp))     # the energy in J/mol
                    
                    Stepis.append(StepID[imtd])
                    Stepts.append(StepID[ts])                       
                    Steptemp.append(temp)                       
                    Steptof.append(np.exp(Espan1*96485/R/temp))     # the energy in J/mol
                    StepE.append(Espan)
                    StepEis.append(Gis)
                    StepEts.append(TS)
                    
                    print('==================================')
                    print('TS=',TS)
                    print('Gis=',Gis)
                    print('ts=',ts)
                    print('imtd=',imtd)
                    print('Espan=', Espan)
                else:
                    print('outside the provided range')            
            else: 
                print('before TS point')
                if imtd < (Glen): 
                    Gis = G_g[imtd]
                    Espan0 = TS-Gis
                    #Espan1 = TS-Gis-abs(Grxn)   # the energy in eV
                    Espan1 = TS-Gis+Grxn   # the energy in eV
                    Espan = Espan0
                    Espanlist.append(Espan)
                    EspanISlist.append(imtd)
                    EspanTSlist.append(ts)
                    StepIS.append(StepID[imtd])
                    StepTS.append(StepID[ts])  
                    RR.append(np.exp(Espan1*96485/R/temp))   # the energy in J/mol
                    
                    Stepis.append(StepID[imtd])
                    Stepts.append(StepID[ts])                       
                    Steptemp.append(temp)                       
                    Steptof.append(np.exp(Espan1*96485/R/temp))   # the energy in J/mol
                    StepE.append(Espan)
                    StepEis.append(Gis)
                    StepEts.append(TS)
                    
                    print('==================================')
                    print('TS=',TS)
                    print('Gis=',Gis)
                    print('ts=',ts)
                    print('imtd=',imtd)
                    print('Espan=', Espan)
                else:
                    print('outside the provided range')
        TSTOF.append(RR)
        
    TStof=[]
    i=0
    for i in range(len(TSTOF)):
        TStof.append(sum(TSTOF[i]))     
    
    print('Espanlist=',Espanlist)
    print('EspanTSlist=',EspanTSlist)
    print('EspanISlist=',EspanISlist)
    Emax=max(Espanlist)
    Emax_id=Espanlist.index(Emax)
    Espan = Emax     # the energy in eV

    EspanISid=EspanISlist[Emax_id]                  
    EspanTSid=EspanTSlist[Emax_id]
    EspanISv=G_g[EspanISid]                  
    EspanTSv=G_g[EspanTSid]

    TOFapprox = kB*temp/h * (np.exp(-Espan*96485/R/temp))    # the energy in J/mol
    TOF = kB*temp/h * ((np.exp(-Grxn*96485/R/temp))-1) / sum(TStof)    # the energy in J/mol
    
    print('Emax=',Emax)
    print('Emax_id=',Emax_id)
    print('EspanTSid=',EspanTSid)
    print('EspanISid=',EspanISid)
    print('Espan/1000, EspanISid, EspanTSid, StepID=',Espan/1000, EspanISid, EspanTSid, StepID)

    return Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn, 'not-ISx', StepIS, 'not-TSx', StepTS, TOFapprox, TOF   # to make it kJ/mol
    #return Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn, ISx, StepIS, TSx, StepTS, TOFapprox, TOF, Stepis, StepEis, Stepts, StepEts, StepE, Steptof, Steptemp

#================================================================================================
# SOLUTIONS 2 
#================================================================================================
kB=const.Boltzmann # J/K
N=const.Avogadro # mol^-1
h=const.Planck # J.sam nie-7
R=N*kB # const.gas_constant in J/mol/K    

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

Estep = ['R','TSR','RX','TS2','UX+HX']

pt = [ppt1,ppt2,ppt3,ppt4,ppt5]
titl = [title1,title2,title3,title4,title5]
specie_ID = [engg1,engg2,engg3,engg4,engg5]

cmass = 1e-3 # mg
flowout = 1 # semi-batch is 1 and batch is 0

#rtemp = [530,540,550,560,570] # 550K
rtemp = [550,650,750,850,950] # 550K
rtime = 10**(5) #for 5%/20%

ranP = 1 #len(titl)
ranT = len(rtemp)
i=0

for i in range(ranP):
    ppt=pt[i]
    title=titl[i]
    engg=specie_ID[i]
    j=0
    AA=[]; BB=[]; CC=[]; DD=[]; EE=[]; FF=[]; GG=[]; TS=[]; IS=[]; RR=[]; 
    tofap=[]; tof=[]; stepts=[]; stepis=[]; isx=[]; tsx=[]; effspan=[]
     
    for j in range(ranT):
        
        TEMP = rtemp[j]
    
        #Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn = espan_calc(ppt,title,engg,rtime,TEMP,cmass,flowout)
        Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn, ISx, StepIS, TSx, StepTS, TOFapprox, TOF = xtof_calc(ppt,title,engg,rtime,TEMP,cmass,flowout)
        #Espan, EspanISid, EspanISv, EspanTSid, EspanTSv, StepID, Grxn, ISx, StepIS, TSx, StepTS, TOFapprox, TOF, Stepis, StepEis, Stepts, StepEts, StepE, Steptof, Steptemp  = xtof_calc(ppt,title,engg,rtime,TEMP,cmass,flowout)
        AA.append(Espan*96.485)    # the energy in kJ/mol
        BB.append(EspanISid)
        CC.append(EspanTSid)
        DD.append(TEMP)
        EE.append(StepID[EspanISid])
        FF.append(StepID[EspanTSid])
        GG.append(Grxn*96.485)   # the energy in kJ/mol
        IS.append(EspanISv*96.485)
        TS.append(EspanTSv*96.485)
        RR.append(kB*TEMP/h*np.exp(-Espan*96485/R/TEMP))   # the energy in J/mol
        tofap.append(TOFapprox)
        tof.append(TOF)
        effspan.append(-R*TEMP*np.log(TOF*h/kB/TEMP))
        stepts.append(StepTS)
        stepis.append(StepIS)
        tsx.append(TSx)
        isx.append(ISx)

    Gstd = stdev(GG)

    Astd = stdev(AA)
    Amx = mean(AA)    

    Rstd = stdev(RR)
    Rmx = mean(RR)   

    csvfile = "/home/toyegoke/ENS_kinetics/Final-espan/"+title+"_data_espan2_.csv"  
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['E-span(kJ/mol)', 'TS_step', 'TDTS(kJ/mol)', 'TD-TS', 'IS_step', 'TDIS(kJ/mol)', 'TD-IS', 'Grxn(kJ/mol)', 'Temp in K', 'Rate(TOF)'])
        
        for j in range(ranT):            
            wr.writerow([ AA[j], CC[j], TS[j], FF[j], BB[j], IS[j], EE[j], GG[j], DD[j], RR[j] ])
        wr.writerow(['E-span-Stdev', Astd, 'E-span-Mean', Amx, 'Grxn', Gstd, 'Rate(TOF)', Rstd])
        
        wr.writerow(['Temp_in_K','TOF_approx','Rate(TOF)','Original_TOF','Energy_Span(kJ/mol)','Grxn(kJ/mol)','effspan(kJ/mol)'])
        for j in range(ranT):            
            wr.writerow([ DD[j], tofap[j], RR[j], tof[j], AA[j], GG[j], effspan[j]/1000 ])

#================================================================================================
