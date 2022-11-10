import numpy as np
import matplotlib.pyplot as plt
import glob, os
from astropy.io import fits
import matplotlib.colors as colors
import antlib
from scipy.optimize import curve_fit as cf

def gauss(x, A, sigma, x0):
    return(A * np.exp(- (x - x0) ** 2 / sigma))

def gauss3_all_parms(x, A1, A2, A3, sigma1, sigma2, sigma3, x0, deltax):
    
    return(gauss(x, A1, sigma1, x0 + deltax) +
           gauss(x, A2, sigma2, x0) +
           gauss(x, A3, sigma3, x0 - deltax))

def gauss3(x, parms):
    return(gauss(x, parms[0], parms[3], parms[6] + parms[7]) +
           gauss(x, parms[1], parms[4], parms[6]           ) +
           gauss(x, parms[2], parms[5], parms[6] - parms[7]))



filenames = glob.glob('./images_fits/061653IV/*.fits')
# filenames = glob.glob('./images_fits/061415IV/*.fits')
filenames = sorted(filenames, key=os.path.basename)

with fits.open(filenames[0], memmap = True) as f:
    f.verify('silentfix')
    imgI = f[0].data
    freq = float(f[0].header["FREQUENC"])
    delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]

histogram, bin_edges = np.histogram(imgI, bins = 10000)
bin_edges = bin_edges[0:-1]

arg_zero = int(np.argwhere(histogram == histogram.max()))
histogram[arg_zero] = np.mean([histogram[arg_zero + 1], histogram[arg_zero - 1]]) 

histogram = histogram / histogram.max() * 100

bord_l = -1
bord_r = 1

parms, acc = cf(gauss3_all_parms, bin_edges, histogram, 
                bounds = ([0, 0, 0, 10 ** (-5), 10 ** (-5), 10 ** (-5), bord_l, 0.25],
                [histogram.max() * 1.1, histogram.max() * 1.1, histogram.max() * 1.1,
                10, 10, 10, bord_r, 5]))

# parms, acc = cf(gauss3_all_parms, bin_edges, histogram, 
#                 bounds = ([0, 0, 0, 10 ** (-5), 10 ** (-5), 10 ** (-5), -1, -1, -1],
#                 [histogram.max() * 1.1, histogram.max() * 1.1, histogram.max() * 1.1,
#                 10, 10, 10, 1, 1, 1]))

print(f"A0 = {parms[0]}, A1 = {parms[1]}, A2 = {parms[2]}")
print(f"sig0 = {parms[0]}, sig1 = {parms[1]}, sig2 = {parms[2]}")
print(f"x0 = {parms[6] - parms[7]}, delta x = {parms[7]}")
print(f"x0 = {parms[6]}, delta x = {parms[7]}")

plt.figure()
plt.scatter(bin_edges, histogram)
plt.plot(bin_edges, gauss3(bin_edges, parms), lw = 4, color = "red")
plt.axvline(x = parms[6], linestyle = '--', color = "black")
plt.axvline(x = parms[6] + parms[7], linestyle = '--', color = "black")
plt.axvline(x = parms[6] - parms[7], linestyle = '--', color = "black")
plt.axis([-2, 2, 0, 110])

plt.figure()
plt.imshow(imgI, origin = "lower", norm = colors.Normalize(vmin=-2, vmax=2))
# plt.imshow(imgI, origin = "lower")
plt.show()


