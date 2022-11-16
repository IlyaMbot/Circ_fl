import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as anime
import os, glob
from astropy.io import fits
from scipy.ndimage import gaussian_filter, center_of_mass

'''
def find_centre_smoothly(orig, shape, sigma = 1, cycles = 1):
        """
        Search for "smooth" coordinates of the brightest pixel  
        """
        
        x, y = np.argwhere(orig == orig.max())[0][0], np.argwhere(orig == orig.max())[0][1]
        frame = orig[ x - shape[0] : x + shape[0], y - shape[1] : y + shape[1]]
        frame = gaussian_filter(frame, sigma = sigma)
        
        for i in range(cycles):

                x1, y1 = np.argwhere(frame == frame.max())[0][0], np.argwhere(frame == frame.max())[0][1]
                frame = gaussian_filter(frame, sigma = sigma)
                # print(x,y)
        
        x, y = x - shape[0] + x1, y - shape[0] + y1

        return(x, y)
'''

def find_excluding_max(orig, shape, cycles = 1):
        _orig = np.copy(orig)

        _orig[_orig < 0] = 0

        x, y = np.argwhere(_orig == _orig.max())[0][0], np.argwhere(_orig == _orig.max())[0][1]
        _frame = np.copy(_orig[ x - shape[0] : x + shape[0], y - shape[1] : y + shape[1]])

        for i in range(cycles):
                x1, y1 = np.argwhere(_frame == _frame.max())[0][0], np.argwhere(_frame == _frame.max())[0][1]
                _frame[x1,y1] = -1
        
        _frame *= -1
        _frame[_frame < 0] = 0
        x1, y1 = int(center_of_mass(_frame)[0]), int(center_of_mass(_frame)[1])
        
        x, y = x - shape[0] + x1, y - shape[0] + y1

        return(x, y)


timestep = 1.375
startf = 10
endf = 20

frameshape_1 = np.array([30, 30]) 
frameshape_2 = np.array([15, 15]) 

#-----------------------------------------------------------------------------

filenames = glob.glob(f'./images_fits/centredflare/*.fits')
filenames = sorted(filenames, key=os.path.basename)


images = []
frames = []
i = 0

#Read-I-files-----------------------------------------------------------
for filename in filenames:

        with fits.open(filename, memmap = True) as f:
                f.verify('silentfix')
                original = f[0].data
        
        original = original / np.max(original) * 100
        images.append(original)

fig, ax = plt.subplots(figsize = (8,6))
images_cutted = []

for imgI in images:
        #cut active region from image
        img0 = imgI[:250, :250]

        x, y = np.argwhere(img0 == img0.max())[0][0], np.argwhere(img0 == img0.max())[0][1]        
        
        # x1, y1 = int(center_of_mass(frame)[0]), int(center_of_mass(frame)[1])
        # print(x1,y1)
        frame_wfl = imgI[ x - frameshape_2[0]: x + frameshape_2[0], y - frameshape_2[1] : y + frameshape_2[1]]
        images_cutted.append(frame_wfl)
        
        
        #Plot-figures-----------------------------------------------------------------
        size = 16

        # levelsI = np.array( [imgI1.max()/ 100 * i for i in range(beg, 100, step)] )
        
        title = ax.text(0.5,1.05,f"Active region, SRH I", size = size *1.2, ha="center", transform=ax.transAxes, weight = "bold")
        im = plt.imshow(frame_wfl, cmap = 'hot', origin = "lower", norm = colors.Normalize(vmin = -10, vmax = 100), animated = True)
        frames.append([im, title])
        i += 1

np.save("cutted_images", images_cutted)

plt.colorbar()
plt.tight_layout()

# ani = anime.ArtistAnimation(fig, frames, interval = 100, repeat = True, blit = False)
# ani.save("flare_centered.gif")
