from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, glob
import cv2
from datetime import datetime


#------------------------------------------------------------------

timestep = 1.375
intensity1 = np.array([])
intensity2 = np.array([])
intensity_g = np.array([])
intensity_full = np.array([])
time = np.array([], dtype = "datetime64[ms]")

filenames = glob.glob('./images_fits/061653IV/*_I.fits')
filenames = sorted(filenames, key=os.path.basename)

# 20
    

images = []
secs = ["16:27", "16:54"]
for j in range(2):
    centre_I  = [241 + j, 60 + j]
    for i in range(0, 20):
        #Read-I-and-V-files-----------------------------------------------------------

        frame = int(filenames[0][-9:-7])
        date = f'2022_05_22_T06_{filenames[0][-16:-10]}'

        with fits.open(filenames[i + 20 * j], memmap = True) as f:
                f.verify('silentfix')
                imgI = f[0].data
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                
        #cut-a-frame-from-image--------------------------------------------------------------------------

        # coords_pix = np.array([-6 + centre_I[0], 10 + centre_I[0], -6 + centre_I[1], 10 + centre_I[1]])
        coords_pix = np.array([-60 + centre_I[0], 60 + centre_I[0], -60 + centre_I[1], 60 + centre_I[1]])

        imgI_cut = imgI[int(coords_pix[2]) : int(coords_pix[3]), int(coords_pix[0]) : int(coords_pix[1])]
        imgI_cut[imgI_cut < 0] = 0
        images.append(imgI_cut)
        #1 = circular, 2 = jet

        center_coord1 = (-coords_pix[0] + centre_I[0] , -coords_pix[2] + centre_I[1])
        center_coord2 = (-coords_pix[0] + 5 + centre_I[0], -coords_pix[2] + 3 + centre_I[1])
        print(center_coord1, center_coord2)
        axesLength = (3, 1)
        angle = -60

        mask1 = np.zeros_like(imgI_cut)
        mask2 = np.zeros_like(imgI_cut)
        mask1 = cv2.circle(mask1, center_coord1, 3, (255, 255, 255), -1)
        mask2 = cv2.ellipse(mask2, center_coord2, axesLength, angle, 0, 360, (255,255,255), -1)
        mask1 = np.array(mask1) / np.max(mask1)
        mask2 = np.array(mask2) / np.max(mask2)
        print(np.sum(mask2))

        intensity1 = np.append(intensity1, np.sum(mask1*imgI_cut))
        intensity2 = np.append(intensity2, np.sum(mask2*imgI_cut))
        intensity_g = np.append(intensity_g, np.sum(imgI_cut))
        intensity_full = np.append(intensity_full, np.sum(imgI))

        time = np.append(time, np.datetime64(f"2021-05-22T06:{secs[j]}", "ms") +
                                np.timedelta64(f"{int(timestep * 1000 * i )}", "ms")).astype(datetime)

np.save("jet_flux", np.array([time, intensity2]))

#Plot-figures-----------------------------------------------------------------
images = np.array(images)

size = 16

plt.figure()
plt.title(f"Image", size = size * 1.1, weight = "bold")
levelsI = np.arange(50, 100, 20)
plt.imshow(images[22], cmap = 'hot')
plt.contour(images[18], levelsI, colors = [(0, 0, 1 - j/255) for j in range(50, 255, 20)], origin = "lower")

plt.figure()

plt.subplot(121)
plt.title(f"Mask circ. flare", size = size * 1.1, weight = "bold")
plt.imshow(mask1)
plt.title(f"Masks jet", size = size * 1.1, weight = "bold")
plt.subplot(122)
plt.imshow(mask2)

plt.figure(figsize = (8, 4))
plt.title(f"Intensity profile from 06:16:{time[0]}", size = size * 1.1, weight = "bold")
plt.plot(time, intensity1 , label = "Circular flare flux")
plt.plot(time, intensity2, label = "Jet")
plt.plot(time, intensity_g, label = "general")
plt.plot(time, intensity_full, label = "full_sun")
plt.xlabel("time, seconds", size = size)
plt.ylabel("Intensity [arb. val.]", size = size)
plt.legend(fontsize = size * 0.5)
plt.tight_layout()
# plt.savefig(f"./SRH_I_{date}{frame}.png", transparent = False, dpi = 400, bbox_inches = "tight")
plt.show()

