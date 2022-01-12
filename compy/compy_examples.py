# the only real dependency compy has (beyond Python3 itself) is numpy
import numpy as np

# though we will also want to generate some plots:
import matplotlib.pyplot as plt

# import compy:
import _mc_compy

# --- --- --- --- ---

# first we will use the main compy function "mc_compy" to generate and plot a 1-scatter redistribution function:

# the majority of compy functions require a set of Gauss-Legendre quadrature grid-points and weights in order to perform angular integrals,
# fortunately numpy can generate these for us
MU_I_N  = 128 # note this value is set also in mc_compy.h
mu_i, w = np.polynomial.legendre.leggauss(MU_I_N)

# we also need the grid to which the output CDF corresponds; note these constants appear in mc_compy.h as well
N           = 8000 
ETA_BOUND_L = 1.0e-4
ETA_BOUND_U = 1.0e4
eta_grid    = np.logspace(np.log10(ETA_BOUND_L), np.log10(ETA_BOUND_U), N, endpoint = True)
# the "eta" here is the ratio of the final post-scatter (lab/fluid frame) photon energy to the initial, E_f/E_0

T_e = 1.0e10 # the electron temperature, in K
E_0 = 1.0e5  # the pre-scatter photon energy, in eV

# mc_compy returns 3 floats and 1 numpy array
mc_ratio, sa_ratio, sigma, cdf = _mc_compy.mc_compy(1, T_e, E_0, mu_i, w)

# let's deal with the CDF---first, we want to be able to make PDFs out of it, for which we'll use the following functions
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
def cdf_to_pdf(eta_grid, cdf, smooth_num = 20):
    pdf = np.zeros(len(eta_grid)-1)
    for i in range(0, len(pdf)):
        pdf[i] = (cdf[i+1] - cdf[i])/(eta_grid[i+1] - eta_grid[i])
    eta_grid_c = np.zeros(len(pdf))
    for i in range(0, len(pdf)):
        eta_grid_c[i] = 0.5 * (eta_grid[i] + eta_grid[i+1])
    pdf = smooth(pdf, smooth_num)
    pdf = pdf/np.trapz(pdf, eta_grid_c)
    return (eta_grid_c, pdf)
# basically these just difference the CDF and smooth the result---they're good for plotting, though for any application
# where we actually want to compute transition probabilities, it's probably better just to difference the CDF directly;

eta_grid_c, pdf = cdf_to_pdf(eta_grid, cdf)

# let's see what this looks like
plt.plot(eta_grid_c, pdf)
plt.semilogx()
plt.xlim([ETA_BOUND_L, ETA_BOUND_U])
plt.ylim([0.0, 1.1*pdf.max()])
plt.xlabel('E_f/E_0')
plt.show()

# the mean amplification ratio is calculated two ways
print('mc_ratio = ' + repr(mc_ratio)) # statistically
print('sa_ratio = ' + repr(sa_ratio)) # and determinstically
# also, the thermally-averaged total cross section, in units of the Thomson scattering cross section
print('sigma    = ' + repr(sigma))

# NOTE: this section is relegated to an easy-to-turn-off code block because it takes a few minutes to run
if True:
    # compy can also compute n-scatter redistribution functions the same way:
    mc_ratio_1, sa_ratio, sigma, cdf_1 = _mc_compy.mc_compy(1, T_e, E_0, mu_i, w)
    mc_ratio_2, sa_ratio, sigma, cdf_2 = _mc_compy.mc_compy(2, T_e, E_0, mu_i, w)
    mc_ratio_3, sa_ratio, sigma, cdf_3 = _mc_compy.mc_compy(3, T_e, E_0, mu_i, w)

    eta_grid_c, pdf_1 = cdf_to_pdf(eta_grid, cdf_1)
    eta_grid_c, pdf_2 = cdf_to_pdf(eta_grid, cdf_2)
    eta_grid_c, pdf_3 = cdf_to_pdf(eta_grid, cdf_3)

    plt.plot(eta_grid_c, pdf_1, 'r-', label = 'n = 1')
    plt.plot(eta_grid_c, pdf_2, 'g-', label = 'n = 2')
    plt.plot(eta_grid_c, pdf_3, 'b-', label = 'n = 3')
    plt.semilogx()
    plt.xlim([ETA_BOUND_L, ETA_BOUND_U])
    plt.ylim([0.0, 1.1*max([pdf_1.max(), pdf_2.max(), pdf_3.max()])])
    plt.xlabel('E_f/E_0')
    plt.legend(frameon = False, loc = 'best')
    plt.show()

    # the statistically-determined mean amplification ratios take into account n-scattering;
    print('mc_ratio_1 = ' + repr(mc_ratio_1))
    print('mc_ratio_2 = ' + repr(mc_ratio_2))
    print('mc_ratio_3 = ' + repr(mc_ratio_3))

# --- --- --- --- ---

# let's do some example calculations for the semi-analytic mean amplifcation ratios and cross section values

T_e = 1.0e8
E_0 = np.logspace(-1, 7, num = 50, endpoint = True)
sa_ratio = np.zeros(len(E_0))
sa_sigma = np.zeros(len(E_0))
for i in range(0, len(E_0)):
    print(repr(i+1) + '/' + repr(len(E_0)))
    sa_ratio[i] = _mc_compy.sa_ratio(T_e, E_0[i], mu_i, w)
    sa_sigma[i] = _mc_compy.sa_sigma(T_e, E_0[i], mu_i, w)

plt.plot(E_0, sa_ratio)
plt.semilogx()
plt.xlabel('E_0')
plt.ylabel('<E_f/E_0>')
plt.show()

plt.plot(E_0, sa_sigma)
plt.semilogx()
plt.xlabel('E_0')
plt.ylabel('sigma/sigma_T')
plt.show()

T_e = np.logspace(3, 11, num = 50, endpoint = True)
E_0 = 1.0e3
sa_ratio = np.zeros(len(T_e))
sa_sigma = np.zeros(len(T_e))
for i in range(0, len(T_e)):
    print(repr(i+1) + '/' + repr(len(T_e)))
    sa_ratio[i] = _mc_compy.sa_ratio(T_e[i], E_0, mu_i, w)
    sa_sigma[i] = _mc_compy.sa_sigma(T_e[i], E_0, mu_i, w)

plt.plot(T_e, sa_ratio)
plt.loglog()
plt.xlabel('T_e')
plt.ylabel('<E_f/E_0>')
plt.show()

plt.plot(T_e, sa_sigma)
plt.semilogx()
plt.xlabel('T_e')
plt.ylabel('sigma/sigma_T')
plt.show()
