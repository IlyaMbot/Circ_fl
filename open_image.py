from astropy.io import fits
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os, glob

def open_AIA_fits(filename):
    with fits.open(filename, memmap = False) as f:
            f.verify('silentfix')
            img = f[0].data
            size = f[0].header["NAXIS1"]
            delta = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
            centre   = [f[0].header["CRPIX1"], f[0].header["CRPIX2"]]
            t_obs = f[0].header["T_OBS"]
        
    return(img, size, delta, centre, t_obs[:-1])

filename1 = glob.glob(f'./sandbox*/pictures/061653IV/*18_I.fits')
filename2 = glob.glob(f'./sandbox*/pictures/061653IV/*19_I.fits')
filename3 = glob.glob(f'./sandbox*/pictures/061653IV/*0_I.fits')

with fits.open(filename1[0], memmap = True) as f:
                f.verify('silentfix')
                imgI1 = f[0].data 
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I  = [259, 255]

with fits.open(filename2[0], memmap = True) as f:
                f.verify('silentfix')
                imgI2 = f[0].data 
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I  = [259, 255]

with fits.open(filename3[0], memmap = True) as f:
                f.verify('silentfix')
                imgI3 = f[0].data 
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I  = [259, 255]
                
imgI1[imgI1 < 0] = 0
imgI2[imgI2 < 0] = 0
imgI3[imgI3 < 0] = 0

imgI = imgI1 - imgI2 
imgI[ imgI < 0 ] = 0
imgI = imgI / np.max(imgI) * 100

filenames = ["AIA_ref_jp2/aia.lev1.94A_2021-05-22T06 16 47.12Z.image_lev1.fits", 
            "AIA_ref_jp2/aia.lev1.131A_2021-05-22T06 16 42.62Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.171A_2021-05-22T06 16 45.35Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.193A_2021-05-22T06 16 53.96Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.211A_2021-05-22T06 16 45.63Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.304A_2021-05-22T06 16 53.13Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.335A_2021-05-22T06 16 48.62Z.image_lev1.fits",
            "AIA_ref_jp2/aia.lev1.1700A_2021-05-22T06 16 52.72Z.image_lev1.fits"]

for i in range(len(filenames)): 
    name = filenames[i].split('_')
    freq = name[2].split(".")[-1][:-1]

    imgAIA, size, delt_pix_AIA, centre_AIA, t_obs = open_AIA_fits(filenames[i])

    shiftx = - 8
    shifty = 2

    extentI    = [ -centre_I[0] * delt_pix_I[0] + shiftx, centre_I[0] * delt_pix_I[0] + shiftx, -centre_I[1] * delt_pix_I[1] + shifty, centre_I[1] * delt_pix_I[1] + shifty]
    extent_AIA = [ -centre_AIA[0] * delt_pix_AIA[0], (size - centre_AIA[0]) * delt_pix_AIA[0], -centre_AIA[1] * delt_pix_AIA[1] , (size - centre_AIA[1]) * delt_pix_AIA[1] ]

    #Plot-figures-----------------------------------------------------------------

    levelsI = np.arange(50, 100, 20)
    size = 16

    fig, ax = plt.subplots(figsize = (8,6))
    plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
    plt.ylabel("Helioprojective Latitude, [arcsec]", size = size* 1.1 , weight = "bold")
    ax.tick_params(axis='both', labelsize = size)
    plt.imshow(imgAIA, cmap = 'hot', origin = "lower", extent = extent_AIA, norm = colors.PowerNorm(gamma=0.25, vmin = 0, vmax = 8000))
    ax.text(0.5,1.05,f"AIA {freq}$\AA$ {t_obs}", size = size *1.2, ha="center", transform=ax.transAxes, weight = "bold")
    # plt.colorbar()
    
    CS1 = plt.contour(imgI, levelsI, colors = [(0, 0, 1 - j/255) for j in range(50, 255, 20)], origin = "lower", extent = extentI)
    CS2 = plt.contour(imgI3, levelsI, colors = [(0, 1 - j/255, 0) for j in range(50, 255, 20)], origin = "lower", extent = extentI)
    
    # for j in range(len(levelsI)): 
        # CS1.collections[j].set_label(f"I, {levelsI[j]}%" + "$I_{max}$")
        # CS2.collections[j].set_label(f"I, {levelsI[j]}%" + "$I_{max}$")
    
    plt.grid()
    plt.axis([-550, -250, 250, 500])
    # plt.axis([-1000, 1000, -1000, 1000])
    plt.tight_layout()
    # plt.legend(loc = 'upper left')
    plt.savefig(f"./{t_obs}_{freq}_AIA.png", transparent = False, dpi = 400, bbox_inches = "tight")
    