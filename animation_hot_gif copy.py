import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob
from astropy.io import fits
from scipy.ndimage import center_of_mass
from scipy.ndimage import zoom

#constants-------------------------------------------------------------------
size_factor = 10.

frameshape_1 = np.array([30, 30]) * int(size_factor)
frameshape_2 = np.array([15, 15]) * int(size_factor)

#-----------------------------------------------------------------------------

images = np.empty((0 , 130 - 61, 268 - 174), dtype = "float32") # set of images
frames = [] # set of images + titles for animation

filenames = glob.glob(f'./images_fits/centredflare/*.fits')
filenames = sorted(filenames, key=os.path.basename)

#Read-all-files-----------------------------------------------------------
for filename in filenames:

        with fits.open(filename, memmap = True) as f:
                f.verify('silentfix')
                original = f[0].data
        
        #normalize, cut out an active region and fill the set with images
        original = original / np.max(original) * 100
        original = original[61 : 130, 174 : 268]
        # original = original[ 0 : 250, 0 : 250]
        images = np.append(images ,np.array([original]), axis = 0)
        # images.append(zoom(original, size_factor))


#Animation--------------------------------------------------------
fig, ax = plt.subplots(figsize = (8,6))

i = 0
for img0 in images:

        img0 = zoom(img0, size_factor)

        #cut active region from image
        
        x, y = np.argwhere(img0 == img0.max())[0][0], np.argwhere(img0 == img0.max())[0][1]        
        frame_wfl = img0[ x - frameshape_2[0]: x + frameshape_2[0], y - frameshape_2[1] : y + frameshape_2[1]]
                
        #Plot-figures-----------------------------------------------------------------
        size = 16
        
        title = ax.text(0.5, 1.05, f"Active region, SRH I", size = size * 1.2, ha = "center", transform = ax.transAxes, weight = "bold")
        im = plt.imshow(frame_wfl, cmap = 'hot', origin = "lower", norm = colors.Normalize(vmin = -10, vmax = 100), animated = True)
        frames.append([im, title])
        i += 1


plt.colorbar()
plt.tight_layout()

ani = anime.ArtistAnimation(fig, frames, interval = 100, repeat = True, blit = False)
ani.save("flare_centered.gif")

# plt.imshow(images[0],origin = "lower")
# plt.show()