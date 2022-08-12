from astropy.io import fits
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob

timestep = 1.375
startf = 10
endf = 20

#-----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize = (8,6))
images = []

#Read-I-files-----------------------------------------------------------
for i in range(startf, endf):

        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{i}_I.fits')
        filenames = sorted(filenames, key=os.path.basename)

        frame = int(filenames[0][-9:-7])
        date = f'2022_05_22_T06_{filenames[0][-16:-10]}'

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI1 = f[0].data / 10
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I1  = [259, 255]

        '''///'''
        
        filenames = glob.glob(f'./sandbox*/pictures/061653IV/*{frame - 1}_I.fits')

        filenames = sorted(filenames, key=os.path.basename)

        with fits.open(filenames[0], memmap = True) as f:
                f.verify('silentfix')
                imgI2 = f[0].data / 10
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                centre_I2   = [263, 251]

        imgI = imgI1 - imgI2
        centre_I = centre_I1
        shiftx = - 20
        shifty = 7
        extentI = [ -centre_I[0] * delt_pix_I[0], centre_I[0] * delt_pix_I[0], -centre_I[1] * delt_pix_I[1] , centre_I[1] * delt_pix_I[1] ]
        
        #Plot-figures-----------------------------------------------------------------

        beg = 30
        step = 20
        size = 16

        # levelsI = np.array( [imgI1.max()/ 100 * i for i in range(beg, 100, step)] )

        title = ax.text(0.5,1.05,f"SRH I (difference), time 06:16:{27 + i * timestep:.2f} UT", size = size *1.2, ha="center", transform=ax.transAxes, weight = "bold")
        img = plt.imshow(imgI, cmap = 'seismic', origin = "lower", extent = extentI, norm = colors.Normalize(vmin = -1, vmax = 1))
        images.append([img, title])



plt.xlabel("Helioprojective Longtitude, [arcsec]", size = size * 1.1, weight = "bold")
plt.ylabel("Helioprojective Latitude, [arcsec]", size = size* 1.1 , weight = "bold")
ax.tick_params(axis='both', labelsize = size)
plt.grid()
plt.axis([-550, -250, 250, 500])
plt.colorbar()
plt.tight_layout()

ani = anime.ArtistAnimation(fig, images, interval = 1000, repeat = True, blit = False)
ani.save("flare_diff.gif")
