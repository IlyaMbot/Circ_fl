from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, glob
import cv2
from datetime import datetime

#functions---------------------------------------------------------
def find_all(path):
    return(sorted(glob.glob(path), key=os.path.basename))


def read_cleanfits_SRH(filename : str):
    '''
    Open FITS file to get image and parameters.

    Parameters
    ----------
    filename : str
    
    Returns: image, frequency, time of observation
    '''

    with fits.open(filename, memmap = True) as f:
        f.verify('silentfix')
        image = np.array(f[0].data)
        freq = f[0].header["FREQUENC"]
        time = f[0].header["T-OBS"]

    return(image, freq, time)


#------------------------------------------------------------------

timestep = 1.375
time = np.array([], dtype = "datetime64[ms]")

filenames = find_all('./images_fits/centredflare/*.fits')

for filename in filenames:
    image, freq, time = read_cleanfits_SRH(filename)
    print(freq, time)


    #Plot-figures-----------------------------------------------------------------

    size = 16

    plt.figure()
    plt.title(f"Image", size = size * 1.1, weight = "bold")
    levelsI = np.arange(50, 100, 20)
    plt.imshow(image, cmap = 'hot')

    # plt.figure()

    # plt.subplot(121)
    # plt.title(f"Mask circ. flare", size = size * 1.1, weight = "bold")
    # plt.imshow(mask1)
    # plt.title(f"Masks jet", size = size * 1.1, weight = "bold")
    # plt.subplot(122)
    # plt.imshow(mask2)

    # plt.figure(figsize = (8, 4))
    # plt.title(f"Intensity profile from 06:16:{time[0]}", size = size * 1.1, weight = "bold")
    # plt.plot(time, intensity1 , label = "Circular flare flux")
    # plt.plot(time, intensity2, label = "Jet")
    # plt.plot(time, intensity_g, label = "general")
    # plt.plot(time, intensity_full, label = "full_sun")
    # plt.xlabel("time, seconds", size = size)
    # plt.ylabel("Intensity [arb. val.]", size = size)
    # plt.legend(fontsize = size * 0.5)
    # plt.tight_layout()
    # # plt.savefig(f"./SRH_I_{date}{frame}.png", transparent = False, dpi = 400, bbox_inches = "tight")
    plt.show()

