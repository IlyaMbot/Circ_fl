from astropy.io import fits
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob

timestep = 1.375
startf = 17
endf = 20



def open_AIA_fits(filename):
        with fits.open(filename_AIA[0], memmap = False) as f:
                f.verify('silentfix')
                img = f[0].data
                delta = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre   = [f[0].header["CRPIX1"], f[0].header["CRPIX2"]]
                size = f[0].header["NAXIS1"]
        
        return(img, delta, centre, size)

#-----------------------------------------------------------------------------



#Read-I-files-----------------------------------------------------------
for i in range(startf, endf):

        filename_AIA = glob.glob('AIA_ref_jp2/aia.lev1.171A_2021-05-22T06 16 45.35Z.image_lev1.fits')
        imgAIA, delt_pix_AIA, centre_AIA, size = open_AIA_fits(filename_AIA)

        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{i}_I.fits')
        filenames = sorted(filenames, key=os.path.basename)

        frame = int(filenames[0][-9:-7])
        date = f'2022_05_22_T06_{filenames[0][-16:-10]}'

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI1 = f[0].data 
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I1  = [259, 255]

        '''///'''
        
        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{frame - 1}_I.fits')

        filenames = sorted(filenames, key=os.path.basename)

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI2 = f[0].data
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I2   = [263, 251]

        imgI = imgI1 - imgI2
        imgI = imgI / 11.252059396296758 * 250
        centre_I = centre_I1
        shiftx = - 20
        shifty = 7
        extentI = [ -centre_I[0] * delt_pix_I[0], centre_I[0] * delt_pix_I[0], -centre_I[1] * delt_pix_I[1] , centre_I[1] * delt_pix_I[1] ]
        extent_AIA = [ -centre_AIA[0] * delt_pix_AIA[0], (size - centre_AIA[0]) * delt_pix_AIA[0], -centre_AIA[1] * delt_pix_AIA[1] , (size - centre_AIA[1]) * delt_pix_AIA[1] ]

        #Plot-figures-----------------------------------------------------------------
        size = 16
        beg = 20
        step = 20
        levelsI = np.array([40, 60, 80, 100])

        fig, ax = plt.subplots(figsize = (8,6))
        title = ax.text(0.5,1.05,f"AIA 171 $\AA$ & SRH I (diff), T 06:16:{27 + i * timestep:.2f} UT", size = size *1.2, ha="center", transform=ax.transAxes, weight = "bold")
        plt.imshow(imgAIA / np.max(imgAIA) * 100, cmap = 'Greys_r', origin = "lower", extent = extent_AIA, norm = colors.Normalize(vmin=-10, vmax=30))
        # plt.contour(imgI1, levels = levelsI, colors = [(i/255, 0, 0) for i in np.linspace(200, 255, len(levelsI))], origin = "lower", extent = extentI, linestyles = 'dashed', linewidths = 3)
        plt.contour(imgI1, levels = levelsI, colors = [(0, 0, 0) for i in np.linspace(200, 255, len(levelsI))], origin = "lower", extent = extentI, linestyles = 'dashed', linewidths = 3)
        # plt.contour(imgI1, levels = levelsI, colors = [(0, 0, 0) for i in np.linspace(0, 0, len(levelsI))] , origin = "lower", extent = extentI, linestyles = 'dashed', linewidths = 3)
        # plt.contour(imgI, levels = levelsI, cmap = "Wistia", origin = "lower", extent = extentI)  
        plt.contour(imgI, levels = levelsI, colors = [(0, 0, 0) for i in np.linspace(200, 255, len(levelsI))], origin = "lower", extent = extentI)                
        plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
        plt.ylabel("Helioprojective Latitude, [arcsec]", size = size* 1.1 , weight = "bold")
        ax.tick_params(axis='both', labelsize = size)
        # plt.grid()
        plt.axis([-550, -250, 250, 500])
        # plt.colorbar()
        plt.tight_layout()
        plt.savefig(f"./h/SRH_diff_{i}.png", transparent = False, dpi = 400, bbox_inches = "tight")



