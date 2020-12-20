#-*-coding:Utf-8 -*

import numpy as np
from numpy import power as pw
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as cst
from astropy import units as u
from iminuit import Minuit
import scipy.integrate
import data_cov as dc  

########################################
#----USEFUL CONSTANTS AND FUNCTIONS----#
########################################


#Constants
c = cst.c.to(u.Mpc/u.second).value #c en Mpc/s
sigma_B = cst.sigma_sb.value #Stefan Boltzman constant(kg/K^4/s^3)

def fv(Neff):
    return Neff*(7./8)*pw(4./11,4./3)/(1+Neff*(7./8)*pw(4./11,4./3))

def Omega_r(H,T0_cmb,Neff):
    rho_rp = 4*sigma_B/pw(cst.c.value,3)*pw(T0_cmb,4)
    return rho_rp*(1+Neff*(7./8)*pw(4./11,4./3))/(3*pw(H*3.24078e-20,2)/(8.*np.pi*cst.G.value))


#Comoving distance (Mpc)

def E(z,Or,OMat,Ok,ODE,w_0,w_a,H):
    return np.sqrt(OMat*pw(1+z,3)+Or*pw(1+z,4)+Ok*pw(1+z,2)+ODE*pw(1+z,3*(1+w_0+w_a))*np.exp(-3*w_a*z/(1+z)))

def co_dist(z,Or,OMat,Ok,ODE,w_0,w_a,H):
    H_s = H*3.24078e-20 #H (s^-1)
    I=scipy.integrate.quad(lambda x: 1/E(x,Or,OMat,Ok,ODE,w_0,w_a,H),0,z)[0]

    if Ok < 0 : #Ok <0 => k>0
        return c/(H_s*np.sqrt(abs(Ok)))*np.sin(np.sqrt(abs(Ok))*I)
    elif Ok > 0 : #Ok > 0 => k<0
        return c/(H_s*np.sqrt(abs(Ok)))*np.sinh(np.sqrt(abs(Ok))*I)
    else :
        return c/H_s*I

#Luminosity Distance (z;Or,OMat,Ok,ODE,w_0,w_a,H) (Mpc)
def Dist_lum(z,Or,OMat,Ok,ODE,w_0,w_a,H):
   distL = (1+z)*co_dist(z,Or,OMat,Ok,ODE,w_0,w_a,H)
   return distL

def Dist_moduli_model(z,Or,OMat,Ok,ODE,w_0,w_a,H):
    mod = 5*np.log10(Dist_lum(z,Or,OMat,Ok,ODE,w_0,w_a,H))+25
    return mod

#Angular Distance (z;Or,OMat,Ok,ODE,w_0,w_a,H) (Mpc)
def Dist_ang(z,Or,OMat,Ok,ODE,w_0,w_a,H):
   distA = 1/(1+z)*co_dist(z,Or,OMat,Ok,ODE,w_0,w_a,H)
   return distA

#Dv in BAO(z;Or,OMat,Ok,ODE,w_0,w_a,H)
def Dv(z,Or,OMat,Ok,ODE,w_0,w_a,H):
    H_s = H*3.24078e-20 #H (s^-1) 
    H_z = H_s*E(z,Or,OMat,Ok,ODE,w_0,w_a,H) #Hubble parameter H(z)
    Dv_z =pw(pw(1+z,2)*pw(Dist_ang(z,Or,OMat,Ok,ODE,w_0,w_a,H),2)*c*z/H_z,1./3) #Article Formula
    return Dv_z

#Ratio of the baryon to photon momentum density 
def Ratio_bp(z,wb,T0_cmb):
    R_z = 31.5*wb*pw(T0_cmb/2.7,-4)*pw(z/pw(10.,3),-1)
    return R_z
   
def theta_cmb(wb,Or,OMat,Ok,ODE,w_0,w_a,H,T0_cmb,Neff):
    wm = OMat*pw(H/100,2)
    H_s = H*3.24078e-20

    g1 = 0.0783*pw(wb,-0.238)*pw(1+39.5*pw(wb,0.763),-1)
    g2 = 0.560*pw(1+21.1*pw(wb,1.81),-1)
    z_ls = 1048*(1+0.00124*pw(wb,-0.738))*(1+g1*pw(wm,g2))
    R_zls = Ratio_bp(z_ls,wb,T0_cmb)
    
    a_eq = 2.35*pw(10.,-5)*pw(wm,-1)*pw(1-fv(Neff),-1)*pw(T0_cmb/2.7,4)
    z_eq = (1-a_eq)/a_eq
    R_eq = Ratio_bp(z_eq,wb,T0_cmb)
    
    rs = c*2*np.sqrt(3)/3*1/np.sqrt(OMat*pw(H_s,2))*np.sqrt(a_eq/R_eq)*np.log((np.sqrt(1+R_zls)+np.sqrt(R_zls+R_eq))/(1+np.sqrt(R_eq)))
    
    DA = co_dist(z_ls,Or,OMat,Ok,ODE,w_0,w_a,H)

    return rs/DA

def rs_zdrag(wm,wb,T0_cmb):
    z_eq = 2.50*pw(10.,4)*wm*pw(T0_cmb/2.7,-4)
    R_eq = Ratio_bp(z_eq,wb,T0_cmb)

    k_eq = 7.46*pw(10.,-2)*wm*pw(T0_cmb/2.7,-2)

    b1 = 0.313*pw(wm,-0.419)*(1+0.607*pw(wm,0.674))
    b2 = 0.238*pw(wm,0.223)
    z_d = 1291*pw(wm,0.251)*(1+b1*pw(wb,b2))/(1+0.659*pw(wm,0.828))
    R_d = Ratio_bp(z_d,wb,T0_cmb)

    rs = 2/(3*k_eq)*np.sqrt(6/R_eq)*np.log((np.sqrt(1+R_d)+np.sqrt(R_d+R_eq))/(1+np.sqrt(R_eq)))
    
    return rs
    

#############
#----JLA----#
############# 

#CHI2 CALCULATION#
#Calculation of MB
def MB(logM,M_1B,dm):
    Mb = np.array([M_1B]*len(logM))
    dm = np.array([dm]*len(logM))
    test=logM>10
    return Mb+dm*test 

def Dist_moduli_data(data_eta, data_logM, A, M1, dm):
    mod = A @ data_eta - MB(data_logM,M1,dm) 
    return mod

#Chi2 calculation
def chi2_jla(a,b,M1,dm,Or,OMat,Ok,ODE,w_0,w_a,H): 
    H_s = H*3.24078e-20
    A = dc.A_matrix(a,b)  #A = A_0 + a*A_1 - b*A_2
    mu = Dist_moduli_data(dc.eta,dc.LogM,A,M1,dm)
    mu_model= np.array([Dist_moduli_model(z,Or,OMat,Ok,ODE,w_0,w_a,H) for z in dc.z_cmb]) #mu
    C_inv_JLA = np.linalg.inv(dc.mu_cov(A))
    M = mu-mu_model
    chis = M.T @ C_inv_JLA @ M
    return chis


#############
#----CMB----#
#############

#Data
C_cmb = pw(10.,-7)*np.array([[0.79039,-4.0042,0.80608],[-4.0042,66.950,-6.9243],[0.80608,-6.9243,3.9712]]) #Cov matrix from Plank and WMAP
C_inv_cmb = np.linalg.inv(C_cmb)
nu_cmb = np.array([0.022065,0.1199,1.041]) #nu = (Omega_b*h^2, Omega_DM*h^2, 100*Theta_MC)

#Chi2 calculation
def chi2_cmb(wb,Or,OMat,Ok,ODE,w_0,w_a,H,T0_cmb,Neff):
    H_s = H*3.24078e-20 #H (s^-1) 
    wm = OMat*(H/100)**2
    wc = wm-wb

    H_s = 2.2e-18

    theta_ls = theta_cmb(wb,Or,OMat,Ok,ODE,w_0,w_a,H,T0_cmb,Neff)

    nu = np.array([wb,wc,100*theta_ls]) - nu_cmb

    chis = (nu.T).dot(C_inv_cmb).dot(nu)

    return chis


#############
#----BAO----#
#############

#Data
C_inv_bao = np.diag([4444,215156,721487])
d_z_bao = np.array([0.336,0.1126,0.07315]) #dz for z = 0.106, z = 0.35 and z = 0.57 
z_bao = np.array([0.106,0.35,0.57])

#Chi2 calculation
def chi2_bao(wb,Or,OMat,Ok,ODE,w_0,w_a,H):
   H_s = H*3.24078e-20 #H s-1

   wm = OMat*pw(H/100,2) 

   rs_zd = rs_zdrag(wm,wb)

   Dv_z = np.array([Dv(z,Or,OMat,Ok,ODE,w_0,w_a,H) for z in z_bao]) 
   d_z = rs_zd/Dv_z - d_z_bao

   chis = (d_z.T).dot(C_inv_bao).dot(d_z)
   return chis

