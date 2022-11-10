from astropy.io import fits
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import os, glob
import cv2

def remove_out_of_phase(arrays : np.ndarray, threshold : float = 1.0) -> np.ndarray:
    '''
    Remove all values with derivative < threshold.
    Parameters
    ----------
    arr : np.ndarray
    threshold : "float", optional
         The default is 1.0.
    '''
    if len(arrays.shape) == 1:
        arrays = np.array([arrays])

    for arr in arrays:

        diffdata = np.diff(arr)

        outs = np.argwhere(abs(diffdata) >= np.max(abs(diffdata)) * threshold)
        outs = np.reshape(outs , (len(outs)))

        for i in range(outs.shape[0] - 1):
            arr[outs[i] + 1 : outs[i + 1] + 1 ] -= arr[outs[i] + 1] - arr[outs[i]]

        arr[outs[-1] + 1 : ] -= arr[outs[-1] + 1] - arr[outs[-1]]

    return(arrays)

def into_sfu(arr) -> np.ndarray:
    
    #constants
    c = 3 * 10 ** 10            #[cm/s] speed of light
    kb = 1.38 * 10 ** (-16)     # [cgs] Bolzmann constant
    cgs2sfu = 10 ** 19          # SFU/[erg/s/cm^2/Hz]

    quiet_sun1 = -0.40495
    quiet_sun2 = -0.10282
    quiet_sun3 =  0.37754

    freq = 5.600 * 10 ** 9
    temp = 15744 # average brightness temperature for pixel for quiet sun

    
    #normalization
    arr = np.copy(arr)
    arr = temp * ( (arr - quiet_sun1) / (quiet_sun2 - quiet_sun1) + 1) #arr in brightness temperature
    arr -= 2 * temp

    #calculations
    omega_per_pix = (delt_pix_I[0] / 3600.) ** 2
    flux_factor = cgs2sfu
    fi = 2 * kb * freq ** 2 /(c ** 2) * omega_per_pix * flux_factor
    arr *= fi
    return(arr) 

#constants------------------------------------------------------------------

timestep = 1.375
c = 3 * 10 ** 10          #[cm/s] speed of light
kb = 1.38 * 10 ** (-16)  # [cgs] Bolzmann constant
cgs2sfu = 10 ** 19 # SFU/[erg/s/cm^2/Hz]

#jet-data

jet_int = np.array([])
jet_time = np.array([], dtype = "datetime64[ms]")

filenames = glob.glob('./images_fits/061653IV/*.fits')
filenames = sorted(filenames, key=os.path.basename)


images = []
secs = ["16:27", "16:54"]

for j in range(2):
    centre_I  = [241 + j, 60 + j]
    for i in range(0, 20):
        #Read-I-files-----------------------------------------------------------


        with fits.open(filenames[i + 20 * j], memmap = True) as f:
                f.verify('silentfix')
                imgI = f[0].data
                imgI[imgI<0] = 0
                freq = float(f[0].header["FREQUENC"])
                delt_pix_I = [f[0].header["CDELT1"], f[0].header["CDELT2"]]
                
        #cut-a-frame-from-image--------------------------------------------------------------------------

        coords_pix = np.array([-60 + centre_I[0], 60 + centre_I[0], -60 + centre_I[1], 60 + centre_I[1]])

        imgI_cut = imgI[int(coords_pix[2]) : int(coords_pix[3]), int(coords_pix[0]) : int(coords_pix[1])]

        #1 = circular, 2 = jet

        center_coord2 = (-coords_pix[0] + 5 + centre_I[0], -coords_pix[2] + 3 + centre_I[1])
        axesLength = (3, 1)
        angle = -60

        mask2 = np.zeros_like(imgI_cut)
        mask2 = cv2.ellipse(mask2, center_coord2, axesLength, angle, 0, 360, (255,255,255), -1)
        mask2 = np.array(mask2) / np.max(mask2) 

        jet_int = np.append(jet_int, np.sum(mask2 * imgI_cut))
        
        jet_time = np.append(jet_time, np.datetime64(f"2021-05-22T06:{secs[j]}", "ms") +
                                np.timedelta64(f"{int(timestep * 1000 * i )}", "ms")).astype(datetime)


# print(jet_int)
#blank_arrays---------------------------------------------------------------

intensity = np.array([])
time = np.array([], dtype = "datetime64[ms]")
xticks = []

#get-all_filenames-from-dir----------------------------------------------------------

filenames = glob.glob('./images_fits/decay/*.fits')
filenames = sorted(filenames, key=os.path.basename)

for i in range(len(filenames)):

        #get-time-from-current-file----------------------------------------------------

        time = np.append(time, np.datetime64(f"2021-05-22T06:15:08", "ms") +
                                np.timedelta64(f"{int(timestep * 1000 * i )}", "ms")).astype(datetime)
        

        #read-files-and-get-data---------------------------------------------------------
        
        with fits.open(filenames[i], memmap = True) as f:
                f.verify('silentfix')
                img = f[0].data
                pixsize = f[0].header["CDELT1"]
                freq = float(f[0].header["FREQUENC"]) 
        
        # mask1 = np.zeros_like(imgI_cut)
        # mask1 = cv2.circle(mask1, center_coord1, 3, (255, 255, 255), -1)
        # mask1 = np.array(mask1) / np.max(mask1)

        intensity = np.append(intensity, np.sum(img))
           
intensity = remove_out_of_phase(intensity, 0.5)[0]

for i in range(len(time)):
        xticks.append(f"{time[i].minute}:{time[i].second}")

#translation into sfu
flux_jet = into_sfu(jet_int)
flux_flare = into_sfu(intensity)

print(into_sfu(np.zeros(2)))
print(np.max(flux_jet))
print(np.max(flux_flare))


# fi = 2 * kb * freq ** 2 /(c ** 2) * omega_per_pix * flux_factor * 15744 / 1.753333
# print(fi_j)


#Plot-figures-----------------------------------------------------------------
# np.save("jet.npy",np.array([np.array(jet_time, dtype = "datetime64"), jet_int]))
# np.save("flare.npy", np.array([np.array(time, dtype = "datetime64"), intensity]))


size = 16


plt.figure(figsize = (7, 6))
plt.title("Intensity profile 22 May 2021", size = size * 1.1, weight = "bold")
plt.plot(time, flux_flare, label = "Flare", lw = 3)
plt.plot(jet_time, flux_jet, label = f"Jet", lw = 3)
# plt.xticks(time, xticks)
plt.xlabel("Time, [UT]", size = size, weight = "bold")
# plt.ylim(0, 1.05)
plt.tick_params(axis = 'x', rotation = 45, labelsize = size * 0.75)
plt.tick_params(axis='y', labelsize = size )
plt.ylabel("Intensity, [normalized]", size = size, weight = "bold")
plt.legend(fontsize = size * 0.9)
plt.tight_layout()
# plt.savefig(f"./SRH_flux_full_jet.png", transparent = False, dpi = 400, bbox_inches = "tight")
plt.show()