# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:41:33 2020

@author: oyego
"""

#=====================================================================================================
# PROGRAMMER DETAILS
#=====================================================================================================
# Engr. Toyese OYEGOKE
# Laboratoire de Chimie, ENS l'Universite de Lyon, Lyon - France
# ToyeseOyegoke@gmail.com, Toyese.oyegoke@ens-lyon.fr, Toyegoke@abu.edu.ng
#=====================================================================================================
#=====================================================================================================

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import csv



# 'CO' will be link for CONTCAR file
# 'b' will be link for OUTCAR file
# 'rho_r' will be the symmetry of the molecule

#contanst
h=const.Planck # 6.62607015e-34 J Hz^-1 or J.s
kB=const.Boltzmann # 1.380649e-23 J K^-1 or J/K
N=const.Avogadro #6.022140857e+23 mol^-1 or 1/mol
c=const.speed_of_light *1e2 # 29979245800 cm s^-1 or cm/s (convert from m s^-1 to cm s^-1)
p=1.013e5 # Pa
pi=math.pi 
R=N*kB #const.gas_constant # 8.314462618 J mol^-1 K^-1
#pi=const.pi  #mu = 1.660538*10**(-27)# kg/amu

#=====================================================================================================
#=====================================================================================================
# THERMODYNAMICS PROPERTIES OF SPECIES ESTIMATIONS
#=====================================================================================================
#=====================================================================================================

def thermoppt(CONTCAR_linK, OUTCAR_linK, SYMMERY_NO, Temperature_in_K, specie_id, ANY_OTHER_LINK = None):
    """
    STATUS = adsorbedsurfac is 'ads' other-cases are 'nads'
    ANY_OTHER_LINK for a case where Electronic energy is somewhere else
    """
    a1=CONTCAR_linK
    a2=a1
    b1=OUTCAR_linK
    b2=b1
    b3=b1
    b4=b1
    elec=ANY_OTHER_LINK
    T=Temperature_in_K # K
    v=R*T/p


    # DECODES WHETHER A SPECIE IS GAS/ADSORBATE/SURFACE SITE
    #========================================================
    if specie_id == 'X':
        STATUS = 'nads' # catalyst site in form of a specie
        print(specie_id,' is the catalyst site specie')
    else:
        def split(word): 
            return list(word) # defines a splitting function
        yy = split(specie_id)
        if 'X' in yy:
            STATUS = 'ads' # adsorbate specie
            print(specie_id,' is the adsorbate specie')
        else:
            STATUS = 'nads' # gas species
            print(specie_id,' is the gas specie')


    # CONFIG ENTROPY
    #================
    Sconfig = 0#1.39*R


    # EXTRACTION OF DEGREE OF FREEDOM (DOF)
    #=======================================
    file= open(b3, 'r')  # OUTCAR in freq file link
    f= open('column', 'w')
    lines = file.readlines() # in the file link
    i=0
    nop=len(lines)
    dof=[] # for appending the deg of freedom
    # function that filters lines
    for i in range(nop):
        cc=lines[i].split(' ')
        #print('cc(',i,')=',cc)
        wavenumber = 'DOF' 
        if(wavenumber in cc):
            #print('line',i,'=',lines[i])  
            #print('wavenumber',i,'=',(lines[i][29:65]))
            #cuu=[str(lines[i][46:57])] # wavenumber in 1/cm
            #print('cuu (',i,')=',cuu)
            dof.append(float(lines[i][29:65]))
            #print('wv_nm=',wv_nm)
    dof # list of wavenumbers in the file link
    print('dof=',dof)
    

    # PRESENCE OF IMAGINARY FREQUENCY
    #==================================
    file= open(b4, 'r')  # OUTTCAR in freq file link
    f= open('column', 'w')
    lines = file.readlines() # in the file link
    i=0
    nop=len(lines)
    img=[] # for appending the deg of freedom
    # function that filters lines
    for i in range(nop):
        cc=lines[i].split(' ')
        #print('cc(',i,')=',cc)
        wavenumber = 'f/i='
        if(wavenumber in cc):
            #print('line',i,'=',lines[i])  
            #print('wavenumber',i,'=',(lines[i]))
            cuu=[str(lines[i])] # wavenumber in 1/cm
            #print('cuu (',i,')=',cuu)
            img.append(lines[i])
            #print('img=',img)
    img # list of wavenumbers in the file link
    #print('img=',img)
    img_no=len(img)
    print('img_no=',img_no)


    # ELECTRONIC ENERGY link
    #=========================
    i=0
    nop=len(lines)
    img=[] # for appending the deg of free

    if elec is None:
        file= open(b4, 'r')  # OUTTCAR in freq file link
        f= open('column', 'w')
        lines = file.readlines() # in the file link
        i=0
        nop=len(lines)
        elect=[] # for appending the deg of freedom
        # function that filters lines
        for i in range(nop):
            cc=lines[i].split(' ')
            #print('cc(',i,')=',cc)
            wavenumber = 'entropy='
            if(wavenumber in cc):
                #print('line',i,'=',lines[i])  
                #print('umber',i,'=',(lines[i][66:90]))
                elect.append(float(lines[i][66:90]))
                #elect.append(float(lines[i][26:46]))
                #print('wv_nm=',wv_nm)
        elect # list of wavenumbers in the file link in eV
        print ('there is electronic contribution on the surface')
        E_elect=elect[0]*96485 # in J/mol (1ev = 96485 J/mol)
        #print('elect_Eprint(np.exp(0))
    else:
        file= open(b4, 'r')  # OUTTCAR in freq file link
        f= open('column', 'w')
        lines = file.readlines() # in the file link
        i=0
        nop=len(lines)
        elect1=[] # for appending the deg of freedom
        # function that filters lines
        for i in range(nop):
            cc=lines[i].split(' ')
            #print('cc(',i,')=',cc)
            wavenumber = 'entropy='
            if(wavenumber in cc):
                #print('line',i,'=',lines[i])  
                #print('umber',i,'=',(lines[i][66:90]))
                elect1.append(float(lines[i][66:90]))
                #elect.append(float(lines[i][26:46]))
                #print('wv_nm=',wv_nm)
        file= open(elec, 'r')  # OUTTCAR in freq file link
        f= open('column', 'w')
        lines = file.readlines() # in the file link
        i=0
        nop=len(lines)
        # function that filters lines
        for i in range(nop):
            cc=lines[i].split(' ')
            #print('cc(',i,')=',cc)
            wavenumber = 'entropy='
            if(wavenumber in cc):
                #print('line',i,'=',lines[i])  
                #print('umber',i,'=',(lines[i][66:90]))
                elect1.append(float(lines[i][66:90]))
                #elect.append(float(lines[i][26:46]))
                #print('wv_nm=',wv_nm)
        elect=min(elect1) # list of wavenumbers in the file link in eV
        print ('there is electronic contribution on the surface')
        E_elect=elect*96485 # in J/mol (1ev = 96485 J/mol)
        #print('elect_Eprint(np.exp(0))

    

    # TOTAL NO OF ATOMS EXTRACTION
    #================================
    fin = open(a1,'r') # CONTCAR in freq file link
    ff = fin.readlines()
    #print('compound name:', ff[0])
    #print('atoms in the compound:', ff[5], '  freq of the atoms:', ff[6])
    nn=(ff[6]).split()
    nn # no of atoms in the structure distribution
    #print('no of atoms in a structure=',nn)


    # EXTRACTION OF SCALING-FACTOR (SF) FOR ATOM COORDINATES
    #==========================================================
    fin2 = open(a2,'r') # CONTCAR in freq file link
    ff2 = fin2.readlines()
    #print('scaling factor:', ff2[1])
    sf=(ff2[1]).split()
    sf # scaling factor for coordinates
    #print('sf=',sf)


    # EXTRACTION OF CELL X-Y-Z COORDINATES FOR THE ATOMS
    #=======================================================
    #for row 1
    R1=[]
    xc1=ff2[2].split(' ')
    #print(xc1)    # x - cell lattice
    i=0
    noo=len(xc1)
    for i in range(noo):
        xx=len(xc1[i])
        qq=0
        if xx > qq:
            R1.append(float(xc1[i]))
        else:
            pass
    #for row 2
    R2=[]
    yc1=ff2[3].split(' ')
    #print(yc1)    # y - cell lattice
    i=0
    noo=len(yc1)
    for i in range(noo):
        yy=len(yc1[i])
        qq=0
        if yy > qq:
            R2.append(float(yc1[i]))
        else:
            pass
    #for row 3    # z - cell lattice
    R3=[]
    zc1=ff2[4].split(' ')
    #print(zc1)
    i=0
    noo=len(zc1)
    for i in range(noo):
        zz=len(zc1[i])
        qq=0
        if zz > qq:
            R3.append(float(zc1[i]))
        else:
            pass
    #print('R1=',R1,'R2=',R2,'R3=',R3)
    #xc = x-coordinate cell (not scaled-up)
    xc=[]
    xc.append(R1[0])
    xc.append(R2[0])
    xc.append(R3[0])
    #yc = y-coordinate cell (not scaled-up)
    yc=[]
    yc.append(R1[1])
    yc.append(R2[1])
    yc.append(R3[1])
    #zc = z-coordinate cell (not scaled-up)
    zc=[]
    zc.append(R1[2])
    zc.append(R2[2])
    zc.append(R3[2])
    #summary of unscaled value of x,y&z of the cell with scale-up factor
    #print('xc=',xc, 'yc=',yc, 'zc=',zc, 'sf=',sf, 'nn=',nn)


    # COMPUTES THE TOTAL NO OF ATOM
    #===============================
    n=len(nn)
    i=1
    nnf=[]
    for i in range(n):
        nne=int(nn[i])
        nnf.append(nne)
        nnf
        sum1=sum(nnf)
    #print('total no of atom =',sum1)


    # COMPUTES THE EXACT POSITIONS OF THE ATOMS
    #===========================================
    st=sum1 #no of atoms in the structure
    i=1
    i0=2
    dif=3
    ii=i0+(i-1)*dif
    atom_x=[]
    atom_y=[]
    atom_z=[]
    x_c=xc*st
    y_c=yc*st
    z_c=zc*st
    s_f=sf*st
    #print('x_c=',x_c,'y_c=',y_c,'z_c=',z_c,'s_f=',s_f)
    for i in range(st):
        #1condition variables
        check = len('Direct')
        check2= (len(str(ff2[7]))) - 1
        #print('wavenumer=', check)
        #print('check2=',check2)
        #print(ff2[7])
        #2condition variables
        xxc1=str(ff2[i+8])
        null=' \n'
        #check3=len(str(ff2[20]))
        #out3=5
        #print(check3<out3)
        if check is check2:
            # when the CONTCAR carries "DIRECT" coordinate
            #print('xx-coordinate cell:', ff2[i+9])
            xxc=(ff2[i+8]).split(' ')
            #print('xxc=',xxc) # scaling factor for coordinates
            #print('xxc=',xxc) # scaling factor for coordinates
            while '' in xxc: xxc.remove('')  
            while 'T' in xxc: xxc.remove('T')  
            while 'F' in xxc: xxc.remove('F')  
            while 'F\n' in xxc: xxc.remove('F\n')   
            while 'T\n' in xxc: xxc.remove('T\n')   
            #print('xxc_ttttttt=',xxc) # scaling factor for coordinates
            #print(xxc[2])
            x_scaleup=(  float(xxc[0])*float(x_c[ii-2]) + float(xxc[1])*float(x_c[ii-1])  +  float(xxc[2])*float(x_c[ii])  )*float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # x_scaleup=float(xxc[2])
            atom_x.append(x_scaleup)
            #print('yy-coordinate cell:', ff2[i+9])
            yyc=(ff2[i+8]).split(' ')
            yyc # scaling factor for coordinates
            #print('yyc=',yyc) # scaling factor for coordinates
            while '' in yyc: yyc.remove('')  
            while 'T' in yyc: yyc.remove('T')  
            while 'F' in yyc: yyc.remove('F')  
            while 'F\n' in yyc: yyc.remove('F\n')   
            while 'T\n' in yyc: yyc.remove('T\n')   
            #print('yyc_ttttttt=',yyc) # scaling factor for coordinates            #x_scaleup=[i * float(xc*sf) in xxc]
            #x_scaleup=[i * float(xc*sf) in xxc]
            #print(yyc[4])
            y_scaleup=(   float(yyc[0])*float(y_c[ii-2])  +   float(yyc[1])*float(y_c[ii-1])  +  float(yyc[2])*float(y_c[ii])    ) *float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # y_scaleup=float(yyc[4])
            atom_y.append(y_scaleup)
            #print('yy-coordinate cell:', ff2[i+9])
            zzc=(ff2[i+8]).split(' ')
            zzc # scaling factor for coordinates
            #print('zzc=',zzc) # scaling factor for coordinates
            while '' in zzc: zzc.remove('')  
            while 'T' in zzc: zzc.remove('T')  
            while 'F' in zzc: zzc.remove('F')  
            while 'F\n' in zzc: zzc.remove('F\n')   
            while 'T\n' in zzc: zzc.remove('T\n')  
            #print('zzc_ttttttt=',zzc) # scaling factor for coordinates            #x_scaleup=[i * float(xc*sf) in xxc]
            #x_scaleup=[i * float(xc*sf) in xxc]
            #print(zzc[6])
            z_scaleup=(   float(zzc[0])*float(z_c[ii-2]) +    float(zzc[1])*float(z_c[ii-1])  + float(zzc[2])*float(z_c[ii])    )*float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # z_scaleup=float(zzc[6])
            atom_z.append(z_scaleup)
        elif xxc1==null: #or check3<out3:
            StopIteration
        else:
            # when the CONTCAR carries "SELECTIVE DYNAMIC (in a line) DIRECT (in another line)" coordinate
            #print('xx-coordinate cell:', ff2[i+9])
            xxc=(ff2[i+9]).split(' ')
            #print('xxc=',xxc) # scaling factor for coordinates
            while '' in xxc: xxc.remove('')  
            while 'T' in xxc: xxc.remove('T')  
            while 'F' in xxc: xxc.remove('F')  
            while 'F\n' in xxc: xxc.remove('F\n')   
            while 'T\n' in xxc: xxc.remove('T\n')   
            #print('xxc_ttttttt=',xxc) # scaling factor for coordinates
            #x_scaleup=[i * float(xc*sf) in xxc]
            #print(xxc[2])
            x_scaleup=(    float(xxc[0])*float(x_c[ii-2])  +  float(xxc[1])*float(x_c[ii-1])   +   float(xxc[2])*float(x_c[ii])    )*float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # x_scaleup=float(xxc[2])
            atom_x.append(x_scaleup)
            #print('yy-coordinate cell:', ff2[i+9])
            yyc=(ff2[i+9]).split(' ')
            yyc # scaling factor for coordinates
            #print('yyc=',yyc) # scaling factor for coordinates
            while '' in yyc: yyc.remove('')  
            while 'T' in yyc: yyc.remove('T')  
            while 'F' in yyc: yyc.remove('F')  
            while 'F\n' in yyc: yyc.remove('F\n')   
            while 'T\n' in yyc: yyc.remove('T\n')   
            #print('yyc_ttttttt=',yyc) # scaling factor for coordinates            #x_scaleup=[i * float(xc*sf) in xxc]
            #print(yyc[4])
            y_scaleup=(   float(yyc[0])*float(y_c[ii-2])  +   float(yyc[1])*float(y_c[ii-1])  +  float(yyc[2])*float(y_c[ii])    ) *float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # y_scaleup=float(yyc[4])
            atom_y.append(y_scaleup)
            #print('yy-coordinate cell:', ff2[i+9])
            zzc=(ff2[i+9]).split(' ')
            zzc # scaling factor for coordinates
            #print('zzc=',zzc) # scaling factor for coordinates
            while '' in zzc: zzc.remove('')  
            while 'T' in zzc: zzc.remove('T')  
            while 'F' in zzc: zzc.remove('F')  
            while 'F\n' in zzc: zzc.remove('F\n')   
            while 'T\n' in zzc: zzc.remove('T\n')  
            #print('zzc_ttttttt=',zzc) # scaling factor for coordinates            #x_scaleup=[i * float(xc*sf) in xxc]
            #x_scaleup=[i * float(xc*sf) in xxc]
            #print(zzc[6])
            z_scaleup=(   float(zzc[0])*float(z_c[ii-2]) +    float(zzc[1])*float(z_c[ii-1])  + float(zzc[2])*float(z_c[ii])    )*float(s_f[i-1]) * 10**(-10)
            # 10^10 was multiplied to convert angstrom to meter
            # z_scaleup=float(zzc[6])
            atom_z.append(z_scaleup)
    #print('atom_x=',atom_x)   
    #print('atom_y=',atom_y)
    #print('atom_z=',atom_z)


    # MASS OF THE ATOM EXTRACTION
    #=============================
    fin1 = open(b1,'r') # OUTCAR in freq file link
    ff1 = fin1.readlines()
    i=0
    nop=len(ff1)-1
    mm=[] # for appending the masses
    # function that filters wavenumber
    for i in range(nop):
        cc=ff1[i] # calling for each line of the file
        wavenumber = 'Mass of Ions in am'
        if(wavenumber in cc):
            mme=ff1[i+1]
            mme1=mme.split() 
            mme1.remove('POMASS') # removing the POMASS variable
            mme1.remove('=') # removing the equality sign
            mme1 # masses of the atom in the structure
    mm1=mme1 # masses of the atom in the structure 
    nop=len(mm1)
    for i in range(nop):
        uuu=float(mm1[i])
        mm.append(uuu)
    mm # masses of the atom in the structure in g/mol
    #print('mm=',mm)


    # COMPUTES THE TOTAL MASS
    #=========================
    n=len(mm)
    i=0
    prod=[]
    for i in range(n):
        product=float(mm[i])*float(nn[i])
        prod.append(product)
        prod
        #print('prod=',prod)
    tot_m=sum(prod)
    #print('total mass from chemical formula =',tot_m)


    # LIST OF ATOM MASS (DETAIL)
    #============================
    #detail list of atom mass
    i=0
    nop=len(mm)
    #mm_list0=[] # for appending the detail list of atom mass
    mm_list=[] # for appending the detail list of atom mass
    for i in range(nop):
        mm_l=mm[i]
        nn_l=int(nn[i])
        mm_add=([float(mm_l)]*int(nn_l))
        mm_list.extend(mm_add)
    #print(mm_list0[0])
    mm_list
    #print('mm_list=',mm_list)
    tot_mm=sum(mm_list) #summation of the masses in the list
    #print('total mass from list of atoms =',tot_mm)


    #CENTRE OF MASS POSITION
    #========================
    # x-coordinate for COM
    n=len(mm_list)
    i=0
    mm_xx=[]
    for i in range(n):
        productt=float(mm_list[i])*float(atom_x[i])
        mm_xx.append(productt)
        mm_xx
        #print('mm_xx=',mm_xx)
    tot_mm_xx=sum(mm_xx)
    x_com=tot_mm_xx/tot_m
    #print('x_com =',x_com) # in meters
    # y-coordinate for COM
    n=len(mm_list)
    i=0
    mm_yy=[]
    for i in range(n):
        product=float(mm_list[i])*float(atom_y[i])
        mm_yy.append(product)
        mm_yy
        #print('mm_yy=',mm_yy)
    tot_mm_yy=sum(mm_yy)
    y_com=tot_mm_yy/tot_m
    #print('y_com =',y_com) # in meters
    # z-coordinate for COM
    n=len(mm_list)
    i=0
    mm_zz=[]
    for i in range(n):
        product=float(mm_list[i])*float(atom_z[i])
        mm_zz.append(product)
        mm_zz
        #print('mm_zz=',mm_zz)
    tot_mm_zz=sum(mm_zz)
    z_com=tot_mm_zz/tot_m
    #print('z_com =',z_com) # in meters
    #print()
    #print('x_com =',x_com,'   y_com =',y_com,'   z_com =',z_com) # in amstrom (A)


    # CALCULATING MOMENT OF INERTIA TENSORS
    #=======================================
    # for Ixx
    n=len(mm_list)
    i=0
    Ixx=[]
    for i in range(n):
        product=float(mm_list[i])*(((float(atom_y[i])-y_com)**2) + ((float(atom_z[i]-z_com))**2))
        Ixx.append(product)
        Ixx
        #print('Ixx=',Izz)
    tot_Ixx=sum(Ixx)
    I_xx=tot_Ixx
    #print('I_xx =',I_xx)
    # for Iyy
    n=len(mm_list)
    i=0
    Iyy=[]
    for i in range(n):
        product=float(mm_list[i])*(((float(atom_x[i])-x_com)**2) + ((float(atom_z[i]-z_com))**2))
        Iyy.append(product)
        Iyy
        #print('Iyy=',Iyy)
    tot_Iyy=sum(Iyy)
    I_yy=tot_Iyy
    #print('I_yy =',I_yy)
    # for Izz
    n=len(mm_list)
    i=0
    Izz=[]
    for i in range(n):
        product=float(mm_list[i])*(((float(atom_y[i])-y_com)**2) + ((float(atom_x[i]-x_com))**2))
        Izz.append(product)
        Izz
        #print('Izz=',Izz)
    tot_Izz=sum(Izz)
    I_zz=tot_Izz
    #print('I_zz =',I_zz)
    # for Ixy
    n=len(mm_list)
    i=0
    Ixy=[]
    for i in range(n):
        product1=(float(mm_list[i]))*((float(atom_y[i]))-y_com)*((float(atom_x[i]))-x_com)
        Ixy.append(product1)
        Ixy
        #print('Ixy=',Ixy)
    tot_Ixy= (sum(Ixy))
    I_xy=tot_Ixy
    I_yx=I_xy
    #print('I_xy =',I_xy)
    # for Izy
    n=len(mm_list)
    i=0
    Izy=[]
    for i in range(n):
        product1=(float(mm_list[i]))*((float(atom_y[i]))-y_com)*((float(atom_z[i]))-z_com)
        Izy.append(product1)
        Izy
        #print('Izy=',Izy)
    tot_Izy= (sum(Izy))
    I_zy=tot_Izy
    I_yz=I_zy
    #print('I_zy =',I_zy)
    # for Izx
    n=len(mm_list)
    i=0
    Izx=[]
    for i in range(n):
        product1=(float(mm_list[i]))*((float(atom_x[i]))-x_com)*((float(atom_z[i]))-z_com)
        Izx.append(product1)
        Izx
        #print('Izx=',Izx)
    tot_Izx= (sum(Izx))
    I_zx=tot_Izx
    I_xz=I_zx
    #print('I_zx =',I_zx)
    #summary of the moment tensor
    I = [[I_xx, I_yx, I_zx], 
    [I_xy, I_yy, I_zy], 
    [I_xz, I_yz, I_zz]]
    #print('I=',I)


    #CALULATING THE PRINCIPAL MOMENT OF INERTIA
    #===========================================
    #diagonalization of tensors matrix
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.linalg as la
    A = np.array(I)
    #print('A=',A)
    eigvals, eigvecs = la.eig(I)
    #print('eigvals=',eigvals)
    #print('eigvecs=',eigvecs)
    eigvals = eigvals.real
    #print('eigvals=',eigvals)
    #principal of moment of inertia
    I_a=eigvals[0]/N/1000 # values of moment of inertia in kg . 'm^2' unit
    I_b=eigvals[1]/N/1000 # values of moment of inertia in kg . 'm^2' unit
    I_c=eigvals[2]/N/1000 # values of moment of inertia in kg . 'm^2' unit
    print('Ia=',I_a,'  Ib=',I_b,'  Ic=',I_c)
    

    # EXTRACTION OF VIBRATIONAL WAVENUMBERS
    #=======================================
    file= open(b2, 'r')  # OUTTCAR in freq file link
    f= open('column', 'w')
    lines = file.readlines() # in the file link
    i=0
    nop=len(lines)-18
    wv_nm=[] # for appending the wavenumbers
    # function that filters wavenumber
    for i in range(nop):
        cc=lines[i].split(' ')
        #print('cc(',i,')=',cc)
        wavenumber = 'cm-1'
        if(wavenumber in cc):
            #print('line',i,'=',lines[i])  
            #print('wavenumber',ii,'=',(lines[i][46:57]))
            #cuu=[str(lines[i][46:57])] # wavenumber in 1/cm
            #print('cuu (',i,')=',cuu)
            wvv=float(lines[i][46:57])
            if wvv > 100:
                wv_nm.append(float(lines[i][46:57]))
            else:
                wv_nm.append(float(100)) # low freq (<100) were all modified to 100
            #print('wv_nm=',wv_nm)
    wv_nm # list of wavenumbers in the file link
    #print('wv_nm=',wv_nm)
    #selection of wavenumber by structure type
    i=0
    wv=[]

    #for classification into MOLECULES/SURFACE
    
    if sum(dof)==0:
        wv = [0] 
        print ('no vibrational effect on the surface')
    else:
        print ('there is vibrational effect on the surface')
        int1=int(3*sum1)
        int2=int(dof[0])
        print(int1)
        print(int2)
        
        Izero=0.0
        Iaa=I_a*(6.022*10**22)*(10**20) #to au.A^2
        Ibb=I_b*(6.022*10**22)*(10**20) #to au.A^2
        Icc=I_c*(6.022*10**22)*(10**20) #to au.A^2
        #print('iAA=',Iaa)
        rt=1*10*(-8)
        at=1*10*(-8)
        Itest1=np.isclose(Izero,Iaa)
        Itest2=np.isclose(Izero,Ibb)
        Itest3=np.isclose(Izero,Icc)
        #print(Itest1,Itest2,Itest3)
    
        if int1==int2:
            #print('molecules')
            if img_no != 0:
                #print('Minima')
                if Itest1 or Itest2 or Itest3:
                    i=0
                    print('molecules-minima-linear')
                    nop=(len(wv_nm))-5 # i.e. 3N-5
                    for i in range(nop):
                        wv.append(wv_nm[i])
                        #print(wv)
                else:
                    i=0
                    print('molecules-minima-non-linear')
                    nop=(len(wv_nm))-6 # i.e. 3N-6
                    for i in range(nop):
                        wv.append(wv_nm[i])
                        #print(wv)
            else:
                #print('TS')
                if Itest1 or Itest2 or Itest3:
                    if img_no>5:
                        dd=img-5
                        i=0
                        print('molecules-TS-linear-with_img')
                        nop=(len(wv_nm))-5-dd # i.e. 3N-5
                        for i in range(nop):
                            wv.append(wv_nm[i])
                            #print(wv)
                    else:
                        i=0
                        print('molecules-TS-linear-no-img')
                        nop=(len(wv_nm))-5 # i.e. 3N-5
                        for i in range(nop):
                            wv.append(wv_nm[i])
                            #print(wv)
                else:
                    if img_no>6:
                        dd=img-6    
                        i=0
                        print('molecules-TS-nonlinear-with_img')
                        nop=(len(wv_nm))-6-dd # i.e. 3N-6
                        for i in range(nop):
                            wv.append(wv_nm[i])
                            #print(wv)
                    else:
                        i=0
                        print('molecules-TS-nonlinear-no-img')
                        nop=(len(wv_nm))-6 # i.e. 3N-6
                        for i in range(nop):
                            wv.append(wv_nm[i])
                            #print(wv)
        else:
            #print('surfaces')
            if img_no>0:
                dd=img_no    
                print('surfaces-with_img')
                i=0
                nop=int(dof[0])-dd
                for i in range(nop):
                    wv.append(wv_nm[i])
                    #print(wv)
            else:
                print('surfaces-without_img')
                i=0
                nop=int(dof[0])
                for i in range(nop):
                    wv.append(wv_nm[i])
                    #print(wv)

    #print('wv=',wv)

 
    #TRANSLATION MOTION EFFECT
    #==========================
    #inputs
    rho_r=SYMMERY_NO #float(input('rho_r (symmetry number) = '))
    m=(float(tot_mm))/N/1000 #kg

    #for classification into MOLECULES/SURFACE
    if sum(dof)==0:
        s_t=0
        q_t=0
        cp_t=0
        dH_t=0
        print ('no translational effect on the surface')
    else:    
        print ('there is translational effect on the surface')
        int1=int(3*sum1)
        int2=int(dof[0])
        #print(int1)
        #print(int2)
        if int1==int2:
            #print('molecules')
            #for s_t
            s_t=R*((3/2)*(np.log(2*pi*m/h/h)) + (5/2)*(np.log(kB*T)) - (np.log(p)) + 5/2)
            #P=1 #bar
            #s_t=R*( (5/2)*np.log(T)   +   (3/2)*np.log(tot_mm)  - np.log(P) - 1.1516)
            #for q_t2
            q_t=((2*pi*m*kB*T)**(3/2))/h/h/h*v
            #for cp_t
            cp_t=(5/2)*R
            #for dH_t
            dH_t=(5/2)*R*T
        else:
            #print('surfaces')
            s_t=0
            q_t=0
            cp_t=0
            dH_t=0
    #summary
    #print('s_t=',s_t,'J/K/mol','  cp_t=',cp_t, 'J/K/mol','  dH_t=',dH_t,'J/mol','  q_t=',q_t,)


    #ROTATIONAL MOTION EFFECT
    #=========================
    #for classification into MOLECULES/SURFACE
    if sum(dof)==0:
        Ia=0
        Ib=0
        Ic=0
        q_r=0
        s_r=0
        cp_r=0
        dH_r=0
        tetha_rot=0
        print ('no rotational effect on the surface')
    else:   
        print ('there is rotational effect on the surface')
        int1=int(3*sum1)
        int2=int(dof[0])
        #print(int1)
        #print(int2)
        if int1==int2:
            #print('molecules')
            if Itest1 or Itest2 or Itest3:
                #print('linear')
                #q_r for linear as q_rl
                I_rl=float(max(I_a,I_b,I_c))  # this will pick only on I for the linear molecules after using Itest to confirm the present of one I that is while are 2 equal I(s)
                B_l=h/8/pi/pi/I_rl
                q_r=kB*T/rho_r/h/B_l
                #print('q_r=',q_r, '   I_rl=',I_rl)
                #s_r for linear as s_rl
                #s_r=R*np.log(np.exp(8*pi*pi*kB*T*I_rl/rho_r/h/h))   
                s_r=R*np.log(T) + R*np.log(I_rl) - R*np.log(rho_r) + (R*np.log(8*pi*pi*kB/h/h) + R)     
                #s_r=R*((np.log(kB*T/rho_r/h/B_l))+1)
                #cp_r for linear as cp_rl
                cp_r=R
                #dH_r for linear as dH_rl
                dH_r=R*T
                #tetha_rot
                tetha_rot=h*h/8/pi/pi/I_rl/kB
            else:
                #print('non-linear polyatomic')
                #q_r for non-linear as q_rn
                Ia=float(I_a)
                Ib=float(I_b)
                Ic=float(I_c)
                q_r=(8*(pi**2)/rho_r/h/h/h)*((Ia*Ib*Ic)**0.5)*((2*pi*kB*T)**(3/2))
                #s_r for non-linear as s_rn
                A=h/8/pi/pi/Ia
                B=h/8/pi/pi/Ib
                C=h/8/pi/pi/Ic
                A1=(2*kB*T*Ia)**0.5
                B1=(2*kB*T*Ib)**0.5
                C1=(2*kB*T*Ic)**0.5
                #s_r=R*np.log(((pi**0.5)*A1*B1*C1/(h/2/pi)/(h/2/pi)/(h/2/pi))/rho_r) + 3*R/2
                #s_r=R*(((3/2)*np.log(kB*T/h)) - ((1/2)*np.log(A*B*C/pi)) - (np.log(rho_r)) + 3/2)
                s_r=(3/2)*R*np.log(T) + (1/2)*R*np.log(Ia*Ib*Ic) - R*np.log(rho_r) + (R*(np.log((8*pi*pi*(2*pi*kB)**(3/2))/h/h/h)) + (3/2)*R)
                #cp_r for non-linear as cp_rn
                cp_r=(3/2)*R
                #dH_r for non-linear as dH_rn
                dH_r=(3/2)*R*T
                #tetha_rot
                tetha_rotA=h*h/8/pi/pi/Ia/kB            
                tetha_rotB=h*h/8/pi/pi/Ib/kB            
                tetha_rotC=h*h/8/pi/pi/Ic/kB
                tetha_rot=(tetha_rotA*tetha_rotB*tetha_rotC)**(1/3)
        else:
            #print('surfaces')
                Ia=0
                Ib=0
                Ic=0
                q_r=0
                s_r=0
                cp_r=0
                dH_r=0
                tetha_rot=0

        
    #summary
    #print('s_r=',s_r,'J/K/mol','  cp_r=',cp_r,'J/K/mol','  dH_r=',dH_r,'J/mol', '  q_r=',q_r,)


    #VIBRATION MOTION EFFECTS
    #=========================
    #inputs
    vi=wv #'vibrational wavenumber in 1/cm
    vi_freq=[] # array of vib freq
    tetha_v=[] # array of tetha for vib
    q_v1=[] # array of vib partition function
    s_v1=[] # array of entropy contrib via vib
    cp_v1=[] # array of cp for vib
    dH_v1=[] # array of enthalpy for vib
    i=0
    
    if sum(dof)==0:
        Ia=0
        Ib=0
        Ic=0
        q_v=0
        s_v=0
        cp_v=0
        dH_v=0
        #print ('no rotational effect on the surface')
    else:   
        #print ('there is rotational effect on the surface')
        nop=len(vi)
        for i in range(nop):
            #for list of vi
            vi_scaleup=vi[i]*1.0015 #scale up factor for B3LYP DFT frequency calculation
            vi_fq=float(float(c)*vi_scaleup) #converts vibrational wavenumber (1/cm) to wave frequency (1/s)
            vi_freq.append(vi_fq)  
            #for list of tetha_v and q_v 
            tetha_v0=float(h)*vi_freq[i]/float(kB)/float(T)
            tetha_v.append(tetha_v0)
            q_v0=(1-np.exp(float(-tetha_v0)))**(-1)
            q_v1.append(q_v0)
            #for list of s_v
            #s_v0=R*tetha_v0*(np.exp(-tetha_v0))/(1 - (np.exp(-tetha_v0)) ) - R*(np.log(1 - np.exp(-tetha_v0) ))
            #s_v0=  R*tetha_v0 / (  (np.exp(tetha_v0)) - 1 )  -  R*( np.log(  1 -  np.exp(-tetha_v0) ))
            s_v0=  R*tetha_v0 / ( (np.exp(tetha_v0))-1 ) - R*( np.log(1-np.exp(-tetha_v0) ))
            #s_v0= - R * (np.log( 1 - np.exp(-tetha_v0) ) ) + R * tetha_v0 * (np.exp(-tetha_v0)) / ( 1 - (np.exp(-tetha_v0)) )
            #s_v0= - R * (np.log( 1 - np.exp(-tetha_v0) ) ) + R * tetha_v0 * (np.exp(-tetha_v0)) / ( 1 - (np.exp(-tetha_v0)) )
            s_v1.append(s_v0)   
            #for list of cp_v
            cp_v0=R*(tetha_v0**2)*(np.exp(-tetha_v0))/((1-np.exp(-tetha_v0))**2)
            cp_v1.append(cp_v0)
            #for list dH_v (measured enthalpychange from dH(0) to dH(T))
            dH_v0=R*T*(tetha_v0)*(np.exp(-tetha_v0))/(1-np.exp(-tetha_v0))
            dH_v1.append(dH_v0)
            #print('vi_freq=',vi_freq)
            #print('tetha_v=',tetha_v)
            #print(s_v1)
            #print(cp_v1)
            #print(dH_v1)
            #print('theta=', tetha_v)
            #q_v
            q_v=np.prod(q_v1)
            q_v    
            print(s_v1)
            #print('summmBEFORE===',sum(s_v1))
            #print('freq of sv element=',len(s_v1))
            #print('initial sume of wv=',len(wv_nm))
    
            #for s_v
            s_v=sum(s_v1)
            print('========================================>>>>>>>>>>>>>>>>>>>>>>>>')
            print(s_v)
            print('========================================>>>>>>>>>>>>>>>>>>>>>>>>')
    
            #for cp_v
            cp_v=sum(cp_v1)
            #print(cp_v)
    
            #for dH_v
            dH_v=sum(dH_v1)
            #print(dH_v)
        
    
    #summary of ppties
    #print('s_v=',s_v,'J/K/mol','  cp_v=',cp_v,'J/K/mol','  dH_v=',dH_v,'J/mol', '  q_v=',q_v,)
    #print('vi_freq=',vi_freq)
    
    
    # ZPE CALCULATION
    #=================
    wv # list of wavenumbers in cm^-1 
    ZPE0=0.5*sum(wv)# in cm^-1 (with 0.1214 % (i.e. ZPE0*1.001214) correction, it will equal KAMILA own)
    ZPE1=ZPE0*4.5563e-6 # in hartree (1 cm^-1 = 4.558e-6 hartree)
    ZPE2=ZPE0*1.23981e-4 # in eV (1 cm^1 = 27.2114 eV)
    ZPE=ZPE0*11.9627 # in J/mol (1 cm^-1 = 11.9627J/mol or 0.0119627 kJ/mol)
    


    #OVERALL THERMODYNMAIC PROPERTIES
    #=================================
    if sum(dof)==0: 
        s_v=0
        s_r=0
        s_t=0
        dH_v=0
        dH_r=0
        dH_t=0    
        cp_v=0
        cp_r=0
        cp_t=0
        ZPE=0
        Q = q_t*q_r*q_v
        S = s_t + s_r + s_v
        dH = dH_t + dH_r + dH_v
        H = dH + E_elect + ZPE
        Cp = cp_t + cp_r + cp_v
        G =H-T*S
        print ('ONLY ELECTRONIC IS ACCOUNTED FOR SURFACE')
    else:   
        print ('WHEN OTHER CONTRIBUTION ARE ACCOUNTED FOR SURFACE')
        #ZPE
        ZPE # in J/mol
        #entropy(S)
        adsurfac='ads'
        if STATUS==adsurfac:
            S = s_v + Sconfig
        else:
            S = s_t + s_r + s_v
        #enthalpy(dH)
        dH = dH_t + dH_r + dH_v
        H = dH + E_elect + ZPE 
        #specific heat capacity (Cp)
        Cp = cp_t + cp_r + cp_v
        #overall partition function (Q)
        Q = q_t*q_r*q_v
        #gibbs free energy (G)
        #dG = dH - (T * S)
        #G=dG+E_elect
        #G=H-T*S
        G=H-T*S
    

    """
    #display of results
    print('E_elect(eV)=',E_elect/96485) 
    print('ZPE(eV)=',ZPE/96485)
    print('s_v(eV/K)=',s_v/96485)
    print('s_r(eV/K)=',s_r/96485)
    print('s_t(eV/K)=',s_t/96485)
    print('S(eV/K)=',S/96485)
    print('dH_v(eV/K)=',dH_v/96485)
    print('dH_r(eV/K)=',dH_r/96485)
    print('dH_t(eV/K)=',dH_t/96485)
    print('dH_sum(eV/K)=',dH/96485)
    print('H(eV/K)=',H/96485)
    print('T*S(eV/K)=',T*S/96485)
    print ('G(in eV)=',G/96485)
    print('S(J/mol/K)=',S)

    """
    print('S=',S,'J/K/mol','    Cp=',Cp,'J/K/mol','    dH=',dH,'J/mol', '    Q=',Q)    #, '    dG=',dG,'J/mol',)
    print('S=',S,'J/K/mol','    Cp=',Cp,'J/K/mol','    dH=',dH,'J/mol' , '    G=',G/96485,'eV',)
    print()
    print('s_v(eV/K)&',s_v/96485,'&s_r(eV/K)&',s_r/96485,'&s_t(eV/K)&',s_t/96485,'&E_elect(eV)&',E_elect/96485,'&ZPE(eV)&',ZPE/96485)
    print()

    return S, Cp, dH,  Q,  G, tetha_rot, tot_mm, H, s_t, s_v, s_r, ZPE, E_elect, dH_t, dH_r, dH_v, q_t, q_r, q_v


#=====================================================================================================
#=====================================================================================================
# FOR LOCATING SPECIE DETAILS IN ANY PROPERTIES LIST
#=====================================================================================================
#=====================================================================================================
def specie_index(specie_id,listP):
    """
    Specie_id = code letter
    listP = H, S or G
    """
    no= listP.index(specie_id)
    print(specie_id,' has index ',no)
    return no


#=====================================================================================================
#=====================================================================================================
# FOR COMPUTING ADSORPTION STEP TRANSITION STATE THERMODYNAMIC PROPERTIES
#=====================================================================================================
#=====================================================================================================
def ads_TS(X_gas, Xt_gas):
    """
    X_gas can be either S or H depending properties of interest for gas	
    Xt_gas is the gas translation contribution to the concerned properties which must be similar with X_gas
    """
    X_TS = X_gas - (1/3 * Xt_gas)   # lossing 1/3 of its translation effect to give us the TS-state property
    print ("X_TS = ", X_TS)
    return X_TS    

def adsTS_E(gas_specie_id,surf_specie_id,ads_specie_id, gas_mol, surf_mol, ads_mol, Temperature_in_K, G_m, H_m, dH_t_m, S_m, s_t_m,int_id):
    """
    1-gas_specie_id,
    2-surf_specie_id,
    3-ads_specie_id,
    4-gas_molm,
    5-surf_mol, 
    6-ads_mol
    """
    T= Temperature_in_K
    
    H_gasTS = ads_TS(H_m[specie_index(gas_specie_id,int_id)], dH_t_m[specie_index(gas_specie_id,int_id)])
    S_gasTS = ads_TS(S_m[specie_index(gas_specie_id,int_id)], s_t_m[specie_index(gas_specie_id,int_id)])
    G_gasTS = H_gasTS - T * S_gasTS

    H_gas = H_m[specie_index(gas_specie_id,int_id)]*gas_mol
    S_gas = S_m[specie_index(gas_specie_id,int_id)]*gas_mol
    G_gas = G_m[specie_index(gas_specie_id,int_id)]*gas_mol

    H_surf = H_m[specie_index(surf_specie_id,int_id)]*surf_mol
    S_surf = S_m[specie_index(surf_specie_id,int_id)]*surf_mol
    G_surf = G_m[specie_index(surf_specie_id,int_id)]*surf_mol

    H_ads = H_m[specie_index(ads_specie_id,int_id)]*ads_mol
    S_ads = S_m[specie_index(ads_specie_id,int_id)]*ads_mol
    G_ads = G_m[specie_index(ads_specie_id,int_id)]*ads_mol

    #H_adsTS = H_ads - (H_gas + H_surf)
    #S_adsTS = S_ads - (S_gas + S_surf)
    #G_adsTS = G_ads - (G_gas + G_surf)
    #G_adsTS = H_adsTS - T * S_adsTS
    
    #H_ads_TS = max(H_adsTS,H_gasTS,0)
    #S_ads_TS = max(S_adsTS,S_gasTS,0)
    H_ads_TS = max(H_ads,H_gas+H_surf,H_gasTS+H_surf)
    S_ads_TS = max(S_ads,S_gas+S_surf,S_gasTS+S_surf)
    #S_ads_TS = S_gasTS
    #G_ads_TS = max(G_adsTS,G_gasTS,0)
    G_ads_TS = max(G_ads,G_gas+G_surf,G_gasTS+G_surf)
    #G_ads_TS = min(G_adsTS,G_gasTS,0)
    #G_ads_TS = G_adsTS
    #G_ads_TS = H_ads_TS - T * S_ads_TS
    print('G_ads_TS=',G_ads_TS/96485)
    print('G_ads,G_gas+G_surf,G_gasTS=', G_ads/96485,(G_gas+G_surf)/96485,G_gasTS/96485)
    print('3*G_surf+G_gasTS=', (3*G_surf + G_ads_TS)/96485 )
    print('3*G_surf+G_gas=', (3*G_surf + G_gas)/96485 )
    print('(3*G_surf+G_gasTS)-(3*G_surf + G_gas)=', ((3*G_surf + G_ads_TS)-(3*G_surf + G_gas))/96485 )
    H_m.append(H_ads_TS)
    S_m.append(S_ads_TS)
    G_m.append(G_ads_TS)
    

    
    print('adsorption TS properties according to Campbell et al (2016, ACS) report')
    print('G_gas=', G_gas/96485, 'G_surf=', G_surf/96485, 'G_ads=', G_ads/96485)
    print('G_gas+G_surf=', (G_gas+G_surf)/96485, 'G_ads=', G_ads/96485)
    #print('R+3X=', (G_gas+3*G_surf)/96485)
    #print('TSR+3X=', (G_ads_TS+3*G_surf)/96485)
    #print('RX+2X=', (G_ads+2*G_surf)/96485)
    #print('(TSR+3X)-(R+3X)=', ((G_ads_TS+3*G_surf)/96485) - ((G_gas+3*G_surf)/96485) )
    #print('(RX+2X)-(R+3X)=', ((G_ads+2*G_surf)/96485) - ((G_gas+3*G_surf)/96485) )
    #print('G_gasTS=', G_gasTS/96485, 'G_adsTS=', G_adsTS/96485)
    print('ads_TS=',G_ads_TS/96485)
    #print('H_ads_TS=',H_ads_TS,'S_ads_TS=',S_ads_TS,'G_ads_TS=',G_ads_TS)
    print()    
    
    return

#=====================================================================================================
#=====================================================================================================
# POTENTIAL ENERGY SURFACE SKETCHING FOR DEHYDROGENATION OF 'R' INTO 'P' & H2
#=====================================================================================================
#=====================================================================================================

#WITH ADS/DES TS
def oxacyclopentanol1(energylist):
    
    #Free Gibb (G) energy for the species
    X=energylist[specie_index('X',int_id)]
    R=energylist[specie_index('R',int_id)]
    RX=energylist[specie_index('RX',int_id)]
    UX=energylist[specie_index('UX',int_id)]
    HX=energylist[specie_index('HX',int_id)]
    VX=energylist[specie_index('VX',int_id)]
    PX=energylist[specie_index('PX',int_id)]
    H2=energylist[specie_index('H2',int_id)]
    P=energylist[specie_index('P',int_id)]
    TS1=energylist[specie_index('TS1',int_id)]
    TS2=energylist[specie_index('TS2',int_id)]
    TSR=energylist[specie_index('TSR',int_id)]
    TSP=energylist[specie_index('TSP',int_id)]
    TSH2=energylist[specie_index('TSH2',int_id)]   
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g   


#WITH ADS/DES TS
def oxacyclopentanol2(energylist):
    #Free Gibb (G) energy for the species
    X=energylist[specie_index('X',int_id)]
    R=energylist[specie_index('Rc',int_id)]
    RX=energylist[specie_index('RcX',int_id)]
    UX=energylist[specie_index('UcX',int_id)]
    HX=energylist[specie_index('HX',int_id)]
    VX=energylist[specie_index('VcX',int_id)]
    PX=energylist[specie_index('PcX',int_id)]
    H2=energylist[specie_index('H2',int_id)]
    P=energylist[specie_index('Pc',int_id)]
    TS1=energylist[specie_index('TS1c',int_id)]
    TS2=energylist[specie_index('TS2c',int_id)]
    TSR=energylist[specie_index('TSRc',int_id)]
    TSP=energylist[specie_index('TSPc',int_id)]
    TSH2=energylist[specie_index('TSH2',int_id)]     
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g



#WITH ADS/DES TS
def cyclopentanol(energylist):
    #Free Gibb (G) energy for the species
    X=energylist[specie_index('X',int_id)]
    R=energylist[specie_index('Rb',int_id)]
    RX=energylist[specie_index('RbX',int_id)]
    UX=energylist[specie_index('UbX',int_id)]
    HX=energylist[specie_index('HX',int_id)]
    VX=energylist[specie_index('VbX',int_id)]
    PX=energylist[specie_index('PbX',int_id)]
    H2=energylist[specie_index('H2',int_id)]
    P=energylist[specie_index('Pb',int_id)]
    TS1=energylist[specie_index('TS1b',int_id)]
    TS2=energylist[specie_index('TS2b',int_id)]
    TSR=energylist[specie_index('TSRb',int_id)]
    TSP=energylist[specie_index('TSPb',int_id)]
    TSH2=energylist[specie_index('TSH2',int_id)]     
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g

    

#WITH ADS/DES TS
def propanol(energylist):
    #Free Gibb (G) energy for the species
    X=energylist[specie_index('X',int_id)]
    R=energylist[specie_index('Ra',int_id)]
    RX=energylist[specie_index('RaX',int_id)]
    UX=energylist[specie_index('UaX',int_id)]
    HX=energylist[specie_index('HX',int_id)]
    VX=energylist[specie_index('VaX',int_id)]
    PX=energylist[specie_index('PaX',int_id)]
    H2=energylist[specie_index('H2',int_id)]
    P=energylist[specie_index('Pa',int_id)]
    TS1=energylist[specie_index('TS1a',int_id)]
    TS2=energylist[specie_index('TS2a',int_id)]
    TSR=energylist[specie_index('TSRa',int_id)]
    TSP=energylist[specie_index('TSPa',int_id)]
    TSH2=energylist[specie_index('TSH2',int_id)]     
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g

#WITH ADS/DES TS
def adamantanol1(energylist):
    #Free Gibb (G) energy for the species
    X=energylist[specie_index('X',int_id)]
    R=energylist[specie_index('Rd',int_id)]
    RX=energylist[specie_index('RdX',int_id)]
    UX=energylist[specie_index('UdX',int_id)]
    HX=energylist[specie_index('HX',int_id)]
    VX=energylist[specie_index('VdX',int_id)]
    PX=energylist[specie_index('PdX',int_id)]
    H2=energylist[specie_index('H2',int_id)]
    P=energylist[specie_index('Pd',int_id)]
    TS1=energylist[specie_index('TS1d',int_id)]
    TS2=energylist[specie_index('TS2d',int_id)]
    TSR=energylist[specie_index('TSRd',int_id)]
    TSP=energylist[specie_index('TSPd',int_id)]
    TSH2=energylist[specie_index('TSH2',int_id)]     
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g


def PES(energylist,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2):

    #Free Gibb (G) energy for the species
    X=energylist[specie_index(x,int_id)]
    R=energylist[specie_index(r,int_id)]
    RX=energylist[specie_index(rx,int_id)]
    UX=energylist[specie_index(ux,int_id)]
    HX=energylist[specie_index(hx,int_id)]
    VX=energylist[specie_index(vx,int_id)]
    PX=energylist[specie_index(px,int_id)]
    H2=energylist[specie_index(h2,int_id)]
    P=energylist[specie_index(p,int_id)]
    TS1=energylist[specie_index(ts1,int_id)]
    TS2=energylist[specie_index(ts2,int_id)]
    TSR=energylist[specie_index(tsr,int_id)]
    TSP=energylist[specie_index(tsp,int_id)]
    TSH2=energylist[specie_index(tsh2,int_id)] 
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    
    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    
    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    
    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    
    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 
    
    #TS_transforming_deh2_products
    #ads_prod0=max((VX+2*HX),PX+2*HX)) - (R+3*X)
    #ZP_g.append(ads_prod0/96485)
    #ads_prod0=max((VX+2*HX),PX+2*HX)) - (R+3*X)
    #ZP_g.append(ads_prod0/96485) 
    
    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))
    
    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))
    
    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    return ZP_g
    




def PESdata(T,filename_link_fomat,energylist,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2):

    #Free Gibb (G) energy for the species
    X=energylist[specie_index(x,int_id)]
    R=energylist[specie_index(r,int_id)]
    RX=energylist[specie_index(rx,int_id)]
    UX=energylist[specie_index(ux,int_id)]
    HX=energylist[specie_index(hx,int_id)]
    VX=energylist[specie_index(vx,int_id)]
    PX=energylist[specie_index(px,int_id)]
    H2=energylist[specie_index(h2,int_id)]
    P=energylist[specie_index(p,int_id)]
    TS1=energylist[specie_index(ts1,int_id)]
    TS2=energylist[specie_index(ts2,int_id)]
    TSR=energylist[specie_index(tsr,int_id)]
    TSP=energylist[specie_index(tsp,int_id)]
    TSH2=energylist[specie_index(tsh2,int_id)] 
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)
    gas_phase=(R+3*X) - (R+3*X)
    ZP_g.append(gas_phase/96485)

    #TS_R_ads
    TS_RX=(TSR+2*X) - (R+3*X)  # the TSR already have X present 
    ZP_g.append(TS_RX/96485)
    TS_RX=(TSR+2*X) - (R+3*X)
    ZP_g.append(TS_RX/96485)
    print('(TSR+2*X)= ',int((TSR+2*X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )

    #ads_R
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)
    ads_react=(RX+2*X) - (R+3*X)
    ZP_g.append(ads_react/96485)

    #TS1_deh2
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    TS1_deh2=(TS1+2*X) - (R+3*X)
    ZP_g.append(TS1_deh2/96485)
    
    #deh2_intermediate
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)
    interm=(UX+HX+X) - (R+3*X)
    ZP_g.append(interm/96485)

    #TS2_deh2
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)
    TS2_deh2=(TS2+HX+X) - (R+3*X)
    ZP_g.append(TS2_deh2/96485)

    #deh2_products
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485)
    ads_prod0=(VX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod0/96485) 

    #transformed_ads_product
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    ads_prod=(PX+2*HX) - (R+3*X)
    ZP_g.append(ads_prod/96485)
    print((ads_prod/96485))

    #TS_product_des
    TS_PX=(TSP+2*HX) - (R+3*X) # the TSP already have X present 
    ZP_g.append(TS_PX/96485)
    TS_PX=(TSP+2*HX) - (R+3*X)
    ZP_g.append(TS_PX/96485)
    print('(TSP+2*HX)=', int((TSP+2*HX)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #product_des
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    des_prod=(P+2*HX+X) - (R+3*X)
    ZP_g.append(des_prod/96485)
    print((des_prod/96485))

    #TS_H2_des
    TS_h2=(TSH2+P+X) - (R+3*X)  # the TSH2 already have 2X present 
    ZP_g.append(TS_h2/96485)
    TS_h2=(TSH2+P+X) - (R+3*X)
    ZP_g.append(TS_h2/96485)
    print('(TSH2+P+X)=', int((TSH2+P+X)/96485), '(R+3*X)=',  int((R+3*X)/96485) )
    
    #H2_desorption
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)
    des_h2=(P+H2+3*X) - (R+3*X)
    ZP_g.append(des_h2/96485)

    int_id = ['R+3X','R+3X','TSR+2X','TSR+2X','RX+2X','RX+2X','TS1+2X','TS1+2X','UX+HX+X','UX+HX+X',
              'TS2+HX+X','TS2+HX+X','VX+2HX','VX+2HX','PX+2HX','PX+2HX','TSP+2HX','TSP+2HX','P+2HX+X','P+2HX+X',
              'TSH2+P+X','TSH2+P+X','P+H2+3X','P+H2+3X']

    csvfile = filename_link_fomat 
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(['ID_x'+str(T),str(energylist)+'_x(J/mol)',str(energylist)+'_x(eV)'])
        glen=len(ZP_g)
        for i in range(glen):            
            wr.writerow([int_id[i],ZP_g[i]*96485,ZP_g[i]])

    return ZP_g



def ESPAN_PES(T,filename_link_fomat,energylist,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2):
    
    Gs = PESdata(T,filename_link_fomat,energylist,int_id,x,r,rx,ux,hx,vx,px,h2,p,ts1,ts2,tsr,tsp,tsh2)
    
    #Relative free Gibb (G) energy for the steps
    ZP_g=[]
    
    #gas_phase
    gas_phase=Gs[0]
    ZP_g.append(gas_phase)

    #TS_R_ads
    TS_RX=Gs[2]
    ZP_g.append(TS_RX)

    #ads_R
    ads_react=Gs[4]
    ZP_g.append(ads_react)

    #TS1_deh2
    TS1_deh2=Gs[6]
    ZP_g.append(TS1_deh2)
    
    #deh2_intermediate
    interm=Gs[8]
    ZP_g.append(interm)

    #TS2_deh2
    TS2_deh2=Gs[10]
    ZP_g.append(TS2_deh2)

    #deh2_products
    ads_prod0=Gs[12]
    ZP_g.append(ads_prod0) 
    
    #TS_transforming_deh2_products     JUST ADDED
    ads_prod1=max(Gs[12],Gs[14])
    ZP_g.append(ads_prod1)     

    #transformed_ads_product
    ads_prod=Gs[14]
    ZP_g.append(ads_prod)

    #TS_product_des
    TS_PX=Gs[16]
    ZP_g.append(TS_PX)
    
    #product_des
    des_prod=Gs[18]
    ZP_g.append(des_prod)

    #TS_H2_des
    TS_h2=Gs[20] # the TSH2 already have 2X present 
    ZP_g.append(TS_h2)
    
    #H2_desorption
    des_h2=Gs[22]
    ZP_g.append(des_h2)
    
    step_ID = ['R+3X','TSR+2X','RX+2X','TS1+2X','UX+HX+X','TS2+HX+X','VX+2HX',
              'TSV+2HX','PX+2HX','TSP+2HX','P+2HX+X','TSH2+P+X','P+H2+3X']

    return ZP_g, step_ID
