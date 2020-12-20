# JLA_code

The code need: 

- The data_file and the README : tablef3.dat and README find on cds https://cdsarc.unistra.fr/viz-bin/cat/J/A+A/568/A22
         
- The covmat directory containing C_bias.fits  C_dust.fits  C_model.fits  C_pecvel.fits  sigma_mu.txt
                                                 C_cal.fits   C_host.fits  C_nonia.fits  C_stat.fits
                                                 link : http://supernovae.in2p3.fr/sdss_snls_jla/covmat_v6.tgz 

To use make python3 Fit_jla.py -h (or --help) :

usage: Fit_jla.py [-h] [--alpha ALPHA] [--beta BETA] [--Mb MB] [--dm DM]
                  [--Om OM] [--OL OL] [--Ok OK] [--wb WB] [--w0 W0] [--wa WA]
                  [--H0 H0] [--F_alpha F_ALPHA] [--F_beta F_BETA]
                  [--F_Mb F_MB] [--F_dm F_DM] [--F_Om F_OM] [--F_OL F_OL]
                  [--F_Ok F_OK] [--F_wb F_WB] [--F_w0 F_W0] [--F_wa F_WA]
                  [--F_H0 F_H0] [--U_jla JLA] [--U_cmb CMB] [--U_bao BAO]
                  [--T_cmb T_CMB] [--Neff NEFF] [--U_SDSS U_SDSS]
                  [--U_SNLS U_SNLS] [--U_LOWZ U_LOWZ] [--U_HST U_HST]
                  [--SN_path SN_DATA] [--RM_path README]
                  [--CovDir_path COVDIR] [-P]

Run JLA CHI2 fit

optional arguments:
  -h, --help            show this help message and exit
  -P                    Option for plotting HD (default: False)

SN nuisance parameters:
  --alpha ALPHA         Alpha SALT2 (default: 0.14)
  --beta BETA           Beta SALT2 (default: 3.1)
  --Mb MB               Abs mag SN (default: -19.0)
  --dm DM               Mass step (default: -0.07)

Cosmo parameters:
  --Om OM               OmegaM_0 (default: 0.3)
  --OL OL               OmegaL_0 (default: 0.7)
  --Ok OK               OmegaK_0 (default: 0.0)
  --wb WB               Omega_b_0*h^2 (default: 0.022)
  --w0 W0               rho_L = w0*p (default: -1.0)
  --wa WA               rho_L = (w0+wa(1-a))*p (default: 0.0)
  --H0 H0               Hubble cst (default: 70.0)

Fixed (True)/Free(False) parameters:
  --F_alpha F_ALPHA     Fix Alpha SALT2 (bool) (default: False)
  --F_beta F_BETA       Fix Beta SALT2 (bool) (default: False)
  --F_Mb F_MB           Fix Abs mag SN (bool) (default: False)
  --F_dm F_DM           Fix Mass step (bool) (default: False)
  --F_Om F_OM           Fix OmegaM_0 (bool) (default: False)
  --F_OL F_OL           Fix OmegaL_0 (bool) (default: True)
  --F_Ok F_OK           Fix OmegaK_0 (bool) (default: True)
  --F_wb F_WB           Fix wb (bool) (default: True)
  --F_w0 F_W0           Fix w0 (bool) (default: True)
  --F_wa F_WA           Fix wa (bool) (default: True)
  --F_H0 F_H0           Fix H0 (bool) (default: True)

Chi2 (Used=True):
  --U_jla JLA           Use jla sn (bool) (default: True)
  --U_cmb CMB           Use cmb prior (bool) (default: False)
  --U_bao BAO           Use bao prior (bool) (default: False)

Constants:
  --T_cmb T_CMB         CMB temperature (default: 2.7255)
  --Neff NEFF           Effective nbr of neutrino (default: 3.046)

SN set to use:
  --U_SDSS U_SDSS       Use SDSS set (default: True)
  --U_SNLS U_SNLS       Use SDSS set (default: True)
  --U_LOWZ U_LOWZ       Use LOWZ set (default: True)
  --U_HST U_HST         Use HST set (default: True)

Path to data:
  --SN_path SN_DATA     tablef3.dat path (default: ./data/tablef3.dat)
  --RM_path README      ReadMe path (default: ./data/ReadMe)
  --CovDir_path COVDIR  Cov Matrix dir path (default: ./data/covmat)

