import numpy as np
import matplotlib.pyplot as plt
import glob, os
from astropy.io import fits
import matplotlib.colors as colors
import antlib
from scipy.optimize import curve_fit as cf

def gauss(x, A, sigma, x0):
    return(A * np.exp(- (x - x0) ** 2 / sigma))


freq = 5.600
sfu0 = 15744

filenames = glob.glob('./images_fits/061653IV/*.fits')
# filenames = glob.glob('./images_fits/061415IV/*.fits')
filenames = sorted(filenames, key=os.path.basename)

with fits.open(filenames[0], memmap = True) as f:
    f.verify('silentfix')
    imgI = f[0].data
    freq = float(f[0].header["FREQUENC"])
    delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]

# imgI = imgI[214 : 298, 213 : 298]
images = []
# images.append(imgI)
images.append(imgI[200 : 300, 213 : 300])
images.append(imgI[196 : 315, 299 : 418])
images.append(imgI[344 : 382, 115 : 143])



bins_min, bins_max = 50, 200


for im in images:
    arr_x0 = np.array([])
    for n_bins in range(bins_min, bins_max):

        histogram, bin_edges = np.histogram(im, bins = n_bins)
        bin_edges = bin_edges[0:-1]
        arg_zero = int(np.argwhere(histogram == histogram.max())[0])
        histogram[arg_zero] = np.mean([histogram[arg_zero + 1], histogram[arg_zero - 1]]) 

        histogram = histogram / histogram.max() * 100

        parms, acc = cf(gauss, bin_edges, histogram, bounds = ([0, 10 ** (-5), -1],
                            [100, 10, 1]))

        arr_x0 = np.append(arr_x0, parms[2])

    mean_x0 = arr_x0.mean()
    dispersion = np.mean( [(mean_x0 - x) ** 2 for x in arr_x0 ]  )

    print(f"x0 = {mean_x0:.5f} +- {dispersion:.5f}")

    # print(f"A = {parms[0]} +- {acc[0]}")
    # print(f"sig = {parms[1]} +- {acc[1]}")
    # print(f"x0 = {parms[2]} +- {acc[2]}")

    plt.figure()
    plt.scatter(bin_edges, histogram)
    plt.plot(bin_edges, gauss(bin_edges, parms[0], parms[1], parms[2]), lw = 4, color = "red")
    plt.axvline(x = parms[2], linestyle = '--', color = "black")
    plt.axis([-2, 2, 0, histogram.max() * 1.01])

    plt.figure()
    plt.imshow(im, origin = "lower", norm = colors.Normalize(vmin=-2, vmax=2))
    # plt.imshow(imgI, origin = "lower")
    plt.show()

