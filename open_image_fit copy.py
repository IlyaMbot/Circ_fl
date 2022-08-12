from astropy.io import fits
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob

timestep = 1.375
startf = 10
endf = 20

def turn_into_seconds(hours, minutes, seconds):
        return(hours * 3600 + minutes * 60 + seconds)

def open_AIA_fits(filename):
        with fits.open(filename[0], memmap = False) as f:
                f.verify('silentfix')
                img = f[0].data
                delta = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre   = [f[0].header["CRPIX1"], f[0].header["CRPIX2"]]
        
        return(img, delta, centre)

#-----------------------------------------------------------------------------



#Read-I-and-V-files-----------------------------------------------------------
for i in range(startf, endf):
        # filenames = glob.glob('./sandbox*/pictures/061415/*02_?.fits')
        # filenames = glob.glob('./sandbox*/pictures/061627IV/*19_?.fits')
        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{i}_I.fits')
        filenames = sorted(filenames, key=os.path.basename)

        frame = int(filenames[0][-9:-7])
        date = f'2022_05_22_T06_{filenames[0][-16:-10]}'
        # frame = 20

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI1 = f[0].data / 10
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I1  = [259, 255]

        '''///'''
        
        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{frame - 1}_?.fits')
        # filenames = glob.glob('./sandbox*/pictures/061627IV/*19_?.fits')
        # filenames = glob.glob(f'./sandbox*/pictures/061627IV/*{frame - 1}_?.fits')
        filenames = sorted(filenames, key=os.path.basename)

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI2 = f[0].data / 10
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I2   = [263, 251]

        # imgI1 = np.roll(imgI1, -3  ,axis = 0 )
        # imgI1 = np.roll(imgI1, 2 ,axis = 1 )

        imgI = imgI1 - imgI2
        
        centre_I = centre_I1

        #Read-AIA-files---------------------------------------------------------------
        '''
        mag = True

        if mag == True:
                filename_AIA = glob.glob('./AIA/mag_2021-05-22T06_15_00.fits')
                imgAIA, delt_pix_AIA, centre_AIA = open_AIA_fits(filename_AIA)
                for i in range(2): 
                        imgAIA = np.flip(imgAIA, i)
                imgAIA = np.nan_to_num(imgAIA, nan = 0.0)
                imgAIA[imgAIA > 500] = 500
                imgAIA[imgAIA <- 500] =- 500
                print(imgAIA.max())

        else:
                filename_AIA = glob.glob('./AIA/304_2021-05-22T06_14_05.fits')
                imgAIA, delt_pix_AIA, centre_AIA = open_AIA_fits(filename_AIA)
                imgAIA[imgAIA > 2000] = 2000
                imgAIA[imgAIA < 5] = 5

        '''
        #-----------------------------------------------------------------------------

        shiftx = - 20
        shifty = 7

        extentI    = [ -centre_I[0] * delt_pix_I[0], centre_I[0] * delt_pix_I[0], -centre_I[1] * delt_pix_I[1] , centre_I[1] * delt_pix_I[1] ]
        extentI    = [ -centre_I[0] * delt_pix_I[0], centre_I[0] * delt_pix_I[0], -centre_I[1] * delt_pix_I[1] , centre_I[1] * delt_pix_I[1] ]
        # extent_AIA = [ -centre_AIA[0] * delt_pix_AIA[0], centre_AIA[0] * delt_pix_AIA[0], -centre_AIA[1] * delt_pix_AIA[1] , centre_AIA[1] * delt_pix_AIA[1] ]

        #Plot-figures-----------------------------------------------------------------

        beg = 30
        step = 20

        # t = strftime("%Y_%m_%d_%H:%M:%S", gmtime())
        levelsI = np.array( [imgI1.max()/ 100 * i for i in range(beg, 100, step)] )
        # levelsV = np.array( [imgV.max()/ 100 * i for i in range(beg, 100, 10)] )

        size = 16
        '''
        if mag == True:
                plt.figure(figsize = (20.5,10))
                plt.title(f"Magnetogram HMI", size = size * 1.5)
                plt.imshow(imgAIA, cmap = 'seismic', origin = "lower", extent = extent_AIA)
                plt.colorbar()
                CS1 = plt.contour(imgI, levels = levelsI, colors = [(0, i/255, 0) for i in np.linspace(75, 255, len(levelsI))], origin = "lower", extent = extentI)
                # CS2 = plt.contour(imgV, levels = levelsV, colors = [(0, i/255, 0) for i in np.linspace(70, 255, len(levelsV))], origin = "lower", extent = extentV )
                plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size)
                plt.ylabel("Helioprojective Latitude, [arcsec]", size = size)
                plt.axis([-550, -250, 250, 500])
                # plt.axis([-960, -300, 200, 500])
                plt.grid()
                j = beg
                for i in range(len(levelsI)): 
                        CS1.collections[i].set_label(f"I, {j}%" + "$I_{max}$")
                        j += step
                        # CS2.collections[i].set_label(f"V, {j}%" + "$V_{max}$")

                plt.legend(loc = 'upper left')
                plt.savefig(f"./mag_contour_{t}.png", transparent = False, dpi = 500, bbox_inches = "tight")
        else:
                plt.figure(figsize = (20.5,10))
                plt.title("AIA 304", size = size * 1.5)
                plt.imshow(imgAIA, cmap = 'hot', origin = "lower", extent = extent_AIA, norm = colors.PowerNorm(gamma=0.25))
                plt.colorbar()
                plt.contour(imgI, levels =[87, 97, 107, 207, 307, 407, 507], colors = [(0, 0, i/255) for i in range(50, 255, 20)], origin = "lower", extent = extentI)
                plt.contour(imgV, levels =[87, 97, 107, 207, 307, 407, 507], colors = [(0, i/255, 0) for i in range(50, 255, 20)], origin = "lower", extent = extentV )
                plt.Circle((500, 0), 1000, color='white', fill=True)
                plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size)
                plt.ylabel("Helioprojective Latitude, [arcsec]", size = size)
                plt.axis([-550, -250, 250, 500])
                # plt.axis([-960, -300, 200, 500])
                plt.grid()
                # plt.savefig(f"./AIA_contour_{t}.png", transparent = False, dpi = 500, bbox_inches = "tight")
        '''

        fig, ax = plt.subplots(figsize = (8,6))
        plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
        plt.ylabel("Helioprojective Latitude, [arcsec]", size = size* 1.1 , weight = "bold")
        ax.tick_params(axis='both', labelsize = size)
        ax.text(0.5,1.05,f"SRH I (difference), time 06:16:{27 + i * timestep:.2f} UT", size = size *1.2, ha="center", transform=ax.transAxes, weight = "bold")
        plt.imshow(imgI, cmap = 'seismic', origin = "lower", extent = extentI, norm = colors.Normalize(vmin = -1, vmax = 1))
        # plt.imshow(imgI, cmap = 'hot', origin = "lower", extent = extentI, norm = colors.Normalize(vmin=-500, vmax=8000))
        plt.grid()
        plt.axis([-550, -250, 250, 500])
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(f"./SRH_I_V_{date}{frame}.png", transparent = False, dpi = 400, bbox_inches = "tight")
        



