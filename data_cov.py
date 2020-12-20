#!/usr/bin/python
#Modified (replace pyfits by astropy.io fits)
import numpy as np
import glob
from astropy.io import fits
from astropy.io.ascii import read as aread



def init(fdata,ReadMe,set_use,CovDir):
    #Data
    global data 
    data = aread(fdata,readme=ReadMe)
    selec_set = [data[i]['set'] in set_use for i in range(len(data))]    
    data = data[selec_set]
    global N_SN
    N_SN = len(data)

    #eta vector
    global eta
    eta = np.zeros(3*N_SN) #eta = (mb_1,X1_1,C_1,...,mb_n,X1_n,C_n)
    mb_index=[i%3 == 0 for i in range(3*N_SN)]
    x1_index=[i%3 == 1 for i in range(3*N_SN)]
    c_index=[i%3 == 2 for i in range(3*N_SN)]
    eta[mb_index] = data['mb']
    eta[x1_index] = data['x1'] 
    eta[c_index] = data['c']
    
    #CMB data
    global z_cmb
    global z_hel
    z_cmb, z_hel = data['zcmb'], data['zhel']

    #Stellar mass
    global LogM
    LogM = data['logMst']
    
    #Cov matrix
    global Ceta
    Ceta = sum([fits.open(mat)[0].data for mat in glob.glob(CovDir+'/C*.fits')])
    Ceta = np.array(Ceta)
    selec_cov = [val for val in selec_set for _ in (0, 1, 2)]
    Ceta = Ceta[np.ix_(selec_cov,selec_cov)]

    #M_A matrix
    global M_A
    M_A = np.ones((3,N_SN,3*N_SN))
    i,j = np.ogrid[:N_SN,:3*N_SN]

    for k in range(len(M_A)):
        M_A[k] = M_A[k]*[3*i+k==j] #Def dans l'article mais correction après test

    #Sigma from Eq. 13
    global sigma
    global sigma_pecvel
    sigma = np.loadtxt(CovDir+'/sigma_mu.txt')[selec_set]
    sigma_pecvel = (5 * 150 / 3e5) / (np.log(10.) * sigma[:, 2])

def A_matrix(alpha,beta):
    A = M_A[0] + alpha*M_A[1] - beta*M_A[2]
    return A

def mu_cov(A):
    Cmu = A @ Ceta @ A.T
    # Add diagonal term from Eq. 13
    Cmu[np.diag_indices_from(Cmu)] += sigma[:, 0] ** 2 + sigma[:, 1] ** 2 + sigma_pecvel ** 2 #Rajoute les termes sur la diagonale.
    
    return Cmu
