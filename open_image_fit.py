from astropy.io import fits
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob

timestep = 1.375
startf = 10
endf = 20


def turn_into_seconds(hours, minutes, seconds):
        return(hours * 3600 + minutes * 60 + seconds)

def open_AIA_fits(filename):
        with fits.open(filename_AIA[0], memmap = False) as f:
                f.verify('silentfix')
                img = f[0].data
                delta = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre   = [f[0].header["CRPIX1"], f[0].header["CRPIX2"]]
        
        return(img, delta, centre)

#Read-AIA-files---------------------------------------------------------------

mag = True

if mag == True:
        filename_AIA = glob.glob('./AIA/mag_2021-05-22T06_18_00.fits')
        imgAIA, delt_pix_AIA, centre_AIA = open_AIA_fits(filename_AIA)
        for i in range(2): 
                imgAIA = np.flip(imgAIA, i)
        imgAIA = np.nan_to_num(imgAIA, nan = 0.0)
        print(imgAIA.max())

else:
        filename_AIA = glob.glob('./AIA/304_2021-05-22T06_16_41.fits')
        imgAIA, delt_pix_AIA, centre_AIA = open_AIA_fits(filename_AIA)
        imgAIA[imgAIA > 2000] = 2000
        imgAIA[imgAIA < 5] = 5

#Read-I-and-V-files-----------------------------------------------------------

cadres = [12,18,19]

# for i in range(startf, endf):
for i in cadres:

        # filenames = glob.glob('./sandbox*/pictures/061415/*02_?.fits')
        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{i}_?.fits')
        filenames = sorted(filenames, key=os.path.basename)

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI = f[0].data * 100
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I   = [f[0].header["CRPIX1"], f[0].header["CRPIX2"]]

        with fits.open(filenames[1], memmap = True) as f:
                f.verify('silentfix')
                imgV = f[0].data * 100



        #-----------------------------------------------------------------------------

        shiftx = 0
        shifty = 0

        extentI    = [ -centre_I[0] * delt_pix_I[0] + shiftx , centre_I[0] * delt_pix_I[0] + shiftx, -centre_I[1] * delt_pix_I[1] + shifty , centre_I[1] * delt_pix_I[1] + shifty]
        extentV    = [ -centre_I[0] * delt_pix_I[0] + shiftx, centre_I[0] * delt_pix_I[0] + shiftx, -centre_I[1] * delt_pix_I[1] + shifty, centre_I[1] * delt_pix_I[1] + shifty ]
        extent_AIA = [ -centre_AIA[0] * delt_pix_AIA[0], centre_AIA[0] * delt_pix_AIA[0], -centre_AIA[1] * delt_pix_AIA[1] , centre_AIA[1] * delt_pix_AIA[1] ]

        #Plot-figures-----------------------------------------------------------------

        beg = 20
        step = 20

        levelsI = np.array( [imgI.max()/ 100 * i for i in range(beg, 100, step)] )
        # levelsV = np.array( [imgV.max()/ 100 * i for i in range(beg, 100, 10)] )

        size = 16

        if mag == True:
                fig, ax = plt.subplots(figsize = (16,7.8))
                ax.text(0.5,1.05,f"HMI and SRH I, time 06:16:{27 + i * timestep:.2f} UT", size = size *1.3, weight = "bold", ha="center", transform=ax.transAxes)
                ax.tick_params(axis='both', labelsize = size)
                plt.imshow(imgAIA, cmap = 'seismic', origin = "lower", extent = extent_AIA, animated = True, norm = colors.Normalize(vmin=-500, vmax=500))
                plt.colorbar()
                plt.contour(imgI, levels = levelsI, colors = [(0, i/255, 0) for i in np.linspace(75, 255, len(levelsI))], origin = "lower", extent = extentI)
                plt.axis([-550, -250, 250, 500])
                # plt.axis([-960, -300, 200, 500])
                plt.grid()
                plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
                plt.ylabel("Helioprojective Latitude, [arcsec]", size = size * 1.1, weight = "bold")
                
                plt.savefig(f"./mag_contour_{i}.png", transparent = False, dpi = 330, bbox_inches = "tight")

                # CS2 = plt.contour(imgV, levels = levelsV, colors = [(0, i/255, 0) for i in np.linspace(70, 255, len(levelsV))], origin = "lower", extent = extentV )
                

                '''
                j = beg
                for i in range(len(levelsI)): 
                        CS1.collections[i].set_label(f"I, {j}%" + "$I_{max}$")
                        j += step
                        # CS2.collections[i].set_label(f"V, {j}%" + "$V_{max}$")
                '''

                # plt.legend(loc = 'upper left')
                # plt.savefig(f"./mag_contour_{t}.png", transparent = False, dpi = 500, bbox_inches = "tight")
        else:
                fig, ax = plt.subplots(figsize = (16,7.8))
                ax.text(0.5,1.05,f"AIA 304 $\AA$ and SRH I, time 06:16:{27 + i * timestep:.2f} UTC", size = size *1.3, weight = "bold", ha="center", transform=ax.transAxes)
                ax.tick_params(axis='both', labelsize = size)
                plt.imshow(imgAIA, cmap = 'hot', origin = "lower", extent = extent_AIA, norm = colors.PowerNorm(gamma=0.25))
                plt.colorbar()
                plt.contour(imgI, levels = levelsI, colors = [(0, i/255, 0) for i in np.linspace(75, 255, len(levelsI))], origin = "lower", extent = extentI)
                plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
                plt.ylabel("Helioprojective Latitude, [arcsec]", size = size * 1.1, weight = "bold")
                plt.axis([-550, -250, 250, 500])
                # plt.axis([-960, -300, 200, 500])
                plt.grid()
                plt.savefig(f"./AIA_contour_{i}.png", transparent = False, dpi = 330, bbox_inches = "tight")




        '''
        plt.figure(figsize = (15,7))
        plt.subplot(121)
        plt.title("SRH I")
        plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size)
        plt.ylabel("Helioprojective Latitude, [arcsec]", size = size)
        plt.imshow(imgI, cmap = 'hot', origin = "lower", extent = extentI)
        # plt.imshow(imgI, cmap = 'hot', origin = "lower", extent = extentI, norm = colors.LogNorm(vmin=imgI.min(), vmax=imgI.max()))
        plt.grid()

        plt.subplot(122)
        plt.title("SRH V")
        plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size)
        plt.ylabel("Helioprojective Latitude, [arcsec]", size = size)
        plt.imshow(imgV, cmap = 'hot', origin = "lower", extent = extentI)
        plt.grid()

        plt.tight_layout()
        '''

