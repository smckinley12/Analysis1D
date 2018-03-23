# Scott McKinley
# Code to analyze extracted 1D Mosfire spectra

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussianfit(x,A,mu,sig):
    return A*np.exp(-(x-mu)**2/(2*sig**2))

filename_oned = input("FITS file name, 1D spectra: ")
filename_eps = input("FITS file name, eps data: ")
#filename_oned = "SDF_J_ptg1_v5b_J_OII-NB711_035903_OIII-NB973_073302_1D_00.fits"
#filename_eps = "SDF_J_ptg1_v5b_J_OII-NB711_035903_OIII-NB973_073302_eps.fits"
eval(input("Emission center pixel: "))#cen = 829
eval(input("Curve-fit pixel radius: "))#rpix = 15
eval(input("Noise window: "))#noise = 50
hdu_oned = fits.open(filename_oned)
hdu_eps = fits.open(filename_eps)
data_oned = hdu_oned[0].data
npix = len(data_oned)
header_eps = hdu_eps[0].header
lambda_ref = header_eps['CRVAL1']
lambda_delt = header_eps['CDELT1']
wavelengths = np.linspace(lambda_ref,lambda_ref+lambda_delt*npix,npix)
lambda_fit = wavelengths[cen-rpix:cen+rpix]
data_fit = data_oned[cen-rpix:cen+rpix]
lambda_noise = list(wavelengths[cen-rpix-noise:cen-rpix])+list(wavelengths[cen+rpix:cen+rpix+noise])
data_noise = list(data_oned[cen-rpix-noise:cen-rpix])+list(data_oned[cen+rpix:cen+rpix+noise])
init = [1,wavelengths[cen],(rpix/2)*lambda_delt]
optgauss,covgauss = curve_fit(gaussianfit,lambda_fit,data_fit,p0=init)
plt.plot(lambda_fit,data_fit)
plt.plot(lambda_fit,gaussianfit(lambda_fit,*optgauss))
A,mu,sig = optgauss
linepix = 5*sig/lambda_delt
line_minpix = int(round(((mu-2.5*sig-lambda_ref)/lambda_delt)))
line_maxpix = int(round(((mu+2.5*sig-lambda_ref)/lambda_delt)))
line_counts = np.sum(data_oned[line_minpix:line_maxpix+1])
print("line counts =",line_counts)
std = np.std(data_noise)
threesig = 3*std*np.sqrt(linepix)
print("3 sigma noise =",threesig)


