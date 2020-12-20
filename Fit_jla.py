#-*-coding:Utf-8 -*

import numpy as np
from numpy import power as pw
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy import constants as cst
from astropy import units as u
from iminuit import Minuit
import scipy.integrate
import data_cov as dc  
import Chi2_jla as chi
import pandas as pds
import argparse
from tabulate import tabulate

def Pos_param(v):
    v=np.float(v)
    if v<0:
        raise argparse.ArgumentTypeError('This parameter has to be > 0')
    return v

def Density_param(v):
    v = np.float(v)
    if v < 0 or v > 1:
        raise argparse.ArgumentTypeError('Value has to be between 0 and 1')
    return v

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Run JLA CHI2 fit',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

sn_n_group = parser.add_argument_group('SN nuisance parameters')
sn_n_group.add_argument('--alpha', dest='alpha', help='Alpha SALT2', type=float,default=0.140)
sn_n_group.add_argument('--beta', dest='beta', help='Beta SALT2', type=float,default=3.1)
sn_n_group.add_argument('--Mb', dest='Mb', help='Abs mag SN', type=float, default=-19.)
sn_n_group.add_argument('--dm', dest='dm', help='Mass step', type=float, default=-0.07)

cosmo_group = parser.add_argument_group('Cosmo parameters')
cosmo_group.add_argument('--Om', dest='Om', help='OmegaM_0', type=Density_param, default=0.3)
cosmo_group.add_argument('--OL', dest='OL', help='OmegaL_0', type=Density_param, default=0.7)
cosmo_group.add_argument('--Ok', dest='Ok', help='OmegaK_0', type=Density_param, default=0.0)
cosmo_group.add_argument('--wb', dest='wb', help='Omega_b_0*h^2', type=Density_param, default=0.0220)
cosmo_group.add_argument('--w0', dest='w0', help='rho_L = w0*p', type=float, default=-1.)
cosmo_group.add_argument('--wa', dest='wa', help='rho_L = (w0+wa(1-a))*p', type=float, default=0.)
cosmo_group.add_argument('--H0', dest='H0', help='Hubble cst', type=Pos_param, default=70.0)

fixed_group = parser.add_argument_group('Fixed (True)/Free(False) parameters')
fixed_group.add_argument('--F_alpha', dest='F_alpha', help='Fix Alpha SALT2 (bool)', type=str2bool, default=False)
fixed_group.add_argument('--F_beta',dest='F_beta', help='Fix Beta SALT2 (bool)', type=str2bool, default=False)
fixed_group.add_argument('--F_Mb',dest='F_Mb', help='Fix Abs mag SN (bool)', type=str2bool, default=False)
fixed_group.add_argument('--F_dm',dest='F_dm', help='Fix Mass step (bool)', type=str2bool, default=False)
fixed_group.add_argument('--F_Om', dest='F_Om', help='Fix OmegaM_0 (bool)', type=str2bool, default=False)
fixed_group.add_argument('--F_OL', dest='F_OL', help='Fix OmegaL_0 (bool)', type=str2bool, default=True)
fixed_group.add_argument('--F_Ok', dest='F_Ok', help='Fix OmegaK_0 (bool)', type=str2bool, default=True)
fixed_group.add_argument('--F_wb', dest='F_wb', help='Fix wb (bool)', type=str2bool, default=True)
fixed_group.add_argument('--F_w0', dest='F_w0', help='Fix w0 (bool)', type=str2bool, default=True)
fixed_group.add_argument('--F_wa', dest='F_wa', help='Fix wa (bool)', type=str2bool, default=True)
fixed_group.add_argument('--F_H0', dest='F_H0', help='Fix H0 (bool)', type=str2bool, default=True)

chi2_group = parser.add_argument_group('Chi2 (Used=True)')
chi2_group.add_argument('--U_jla', dest='jla', help='Use jla sn (bool)', type=str2bool, default=True)
chi2_group.add_argument('--U_cmb', dest='cmb', help='Use cmb prior (bool)', type=str2bool, default=False)
chi2_group.add_argument('--U_bao', dest='bao', help='Use bao prior (bool)', type=str2bool, default=False)

cst_group = parser.add_argument_group('Constants')
cst_group.add_argument('--T_cmb', dest='T_cmb', help='CMB temperature', type=Pos_param, default=2.7255)
cst_group.add_argument('--Neff', dest='Neff', help='Effective nbr of neutrino', type=Pos_param, default=3.046)

sn_use_group = parser.add_argument_group('SN set to use')
sn_use_group.add_argument('--U_SDSS', dest='U_SDSS', help='Use SDSS set', type=str2bool, default=True)
sn_use_group.add_argument('--U_SNLS', dest='U_SNLS', help='Use SDSS set', type=str2bool, default=True)
sn_use_group.add_argument('--U_LOWZ', dest='U_LOWZ', help='Use LOWZ set', type=str2bool, default=True)
sn_use_group.add_argument('--U_HST', dest='U_HST', help='Use HST set', type=str2bool, default=True)

data_path_group = parser.add_argument_group('Path to data')
data_path_group.add_argument('--SN_path', dest='sn_data', help='tablef3.dat path', type=str, default='./data/tablef3.dat')
data_path_group.add_argument('--RM_path', dest='readme', help='ReadMe path', type=str, default='./data/ReadMe')
data_path_group.add_argument('--CovDir_path', dest='CovDir', help='Cov Matrix dir path', type=str, default='./data/covmat')

plot =  parser.add_argument('-P', dest='Plot', help='Option for plotting HD', action='store_true')

args = parser.parse_args()

if not args.F_OL and not args.F_Ok:
     raise argparse.ArgumentTypeError("OL and Ok can't varry at the same time")

#Use only some sn
set_use=[]
if args.U_SDSS:
    set_use.append(1)
if args.U_SNLS:
    set_use.append(2)
if args.U_LOWZ:
    set_use.append(3)
if args.U_HST:
    set_use.append(4)

dc.init(args.sn_data,args.readme,set_use,args.CovDir) 

##############
#----CHI2----#
##############

def chi2(x):# x=(a=x[0],b=x[1],M1=x[2],dm=x[3],wb=x[4],OMat=x[5],OL=[6],Ok=x[7],w0=x[8],w_a=x[9],H=x[10])

    Or = chi.Omega_r(x[10],args.T_cmb,args.Neff)
    
    if args.F_OL and args.F_Ok: #OL and Ok fixed
        if x[7] != 0.0: #But Ok modified -> No-FlatLCDM param by Ok
            Ok = x[7] 
            ODE = 1-Or-x[5]-x[7]
        elif x[6] != 0.7: #But OL modified -> No-FlatLCDM param by OL
            ODE = x[6]
            Ok = 1-Or-x[5]-x[6]
        else: #Flat-LCDM
            Ok = 0
            ODE = 1-Or-x[5]
            
    elif args.F_Ok: #
         Ok= 1-Or-x[5]-x[6]
         ODE = x[6]
    elif args.F_OL:
         ODE = 1-Or-x[5]-x[7]
         Ok = x[7]

    if args.jla:
        c2_jla = chi.chi2_jla(x[0],x[1],x[2],x[3],Or,x[5],Ok,ODE,x[8],x[9],x[10]) #chi2_jla(a,b,M1,dm,Or,OMat,Ok,ODE,w_0,w_a,H)
    else:
        c2_jla = 0

    if args.cmb:
        c2_cmb = chi.chi2_cmb(x[4],Or,x[5],Ok,ODE,x[8],x[9],x[10],args.T_cmb,args.Neff) #chi_2_cmb(wb,Or,OMat,Ok,ODE,w_0,w_a,H)
    else :
        c2_cmb = 0

    if args.bao :
        c2_bao = chi.chi2_bao(x[4],Or,x[5],Ok,ODE,x[8],x[9],x[10]) #chi2_bao(wb,Or,OMat,Ok,ODE,w_0,w_a,H)
    else :
        c2_bao = 0

    return c2_jla+c2_bao+c2_cmb


#x=(a=x[0],b=x[1],M1=x[2],dm=x[3],wb=x[4],OMat=x[5],OL=[6],Ok=x[7],w0=x[8],w_a=x[9],H=x[10])
names = ["a", "b", "Mb" , "dm", "wb", "OMat","OL","Ok","w_0","w_a","H0"]
     
x0 = np.array([args.alpha, args.beta, args.Mb, args.dm, args.wb, args.Om, args.OL, args.Ok, args.w0, args.wa, args.H0])
fixed = (args.F_alpha, args.F_beta, args.F_Mb, args.F_dm, args.F_wb, args.F_Om,args.F_OL, args.F_Ok,args.F_w0,args.F_wa, args.F_H0)
limits = ((None,None),(None,None),(None,None),(None,None),(0.,1.),(0.,1.),(0.,1.),(0.,1.),(None,None),(None,None),(0.,None))
chimin = Minuit.from_array_func(chi2,x0,limit=limits,name=names, errordef=1, fix=fixed, error=0.0001)

if args.jla==False:
    chimin.fixed['a']=True
    chimin.fixed['b']=True
    chimin.fixed['Mb']=True
    chimin.fixed['dm']=True

#chimin.migrad()
#chimin.hesse()

names_param = ['alpha','beta','Mb','dm','Or','OMat','OL','Ok','w0','wa','H0','Tcmb','Neff']
results = np.concatenate((chimin.values[:4],[chi.Omega_r(chimin.values['H0'],args.T_cmb,args.Neff)],chimin.values[5:],[args.T_cmb],[args.Neff]),axis=0)


if chimin.fixed['OL'] and chimin.fixed['Ok']:
    results[6]=1-results[5]-results[4]
elif chimin.fixed['OL']:
    results[6] = 1-results[5]-results[7]-results[4]
else:
    results[7] = 1-results[5]-results[6]-results[4]

results_table=pds.DataFrame([names_param,results,np.concatenate((chimin.errors[:],['NoErr','NoErr'])),np.concatenate((chimin.fixed[:],['CST','CST']))],['Param', 'Value', 'Error', 'Fixed'])
print(tabulate(results_table.T,headers='keys',tablefmt='psql'))

##########
#--PLOT--#
##########    
        
if args.Plot:
    a,b,M1,dm,Or,OMat,ODE,Ok,w0,wa,H,Tcmb,Neff = results

    A_fit = dc.A_matrix(a,b)
    C_fit = dc.mu_cov(A_fit)
    mu_fit =  chi.Dist_moduli_data(dc.eta,dc.LogM, A_fit, M1, dm)
    err_fit_mu = np.sqrt(np.diag(C_fit))
    z_plot=np.linspace(np.min(dc.data['zcmb'])-0.001,np.max(dc.data['zcmb']+0.1),1000)
    mu_model =  np.array([chi.Dist_moduli_model(z,Or,OMat,Ok,ODE,w0,wa,H) for z in z_plot])
    mu_comp = np.array([chi.Dist_moduli_model(z,Or,OMat,Ok,ODE,w0,wa,H) for z in dc.data['zcmb']])
    Survey_dic = {1: "SDSS",
                  2: "SNLS",
                  3: "LOW_Z",
                  4: "HST"}
    
    fig=plt.figure()
    plt.rcParams.update({'font.size': 15})
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax0 = plt.subplot(gs[0])
    ax0.set_xscale("log")
    ax0.set_ylabel('$\mu$')
    ax0.tick_params(axis='x',which='both', direction='in')
    ax0.plot(z_plot,mu_model,color='black')
    ax1 = plt.subplot(gs[1],sharex=ax0)
    ax1.set_xlabel('z')
    ax1.set_ylabel('$\mu-\mu_{Model}$')
    for i in set_use:
        selec = dc.data['set']==i
        data_set = dc.data[selec]
        mu_set = mu_fit[selec]
        mu_comp_set=mu_comp[selec]
        err_set = err_fit_mu[selec]    
        ax0.errorbar(data_set['zcmb'],mu_set,yerr=err_set,fmt ='.',elinewidth=0.5,ms=2,label=Survey_dic[i])
        ax1.errorbar(data_set['zcmb'],mu_set-mu_comp_set,yerr=err_set,fmt ='.',elinewidth=0.5,ms=2,label=Survey_dic[i])
    ax0.legend()
    ax1.axhline(color='k')
    ax1.set_yticks(np.arange(-0.4,0.5,step=0.2))
    ax1.grid(ls='--',axis='y')
    ax1.set_ylim(-0.5,0.5)
    plt.subplots_adjust(hspace=.0)
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.show()
