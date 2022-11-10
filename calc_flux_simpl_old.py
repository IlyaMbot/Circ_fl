from astropy.io import fits
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import os, glob

def remove_out_of_phase(arrays : np.ndarray, threshold : float = 1.0) -> np.ndarray:
    '''
    Remove all out of phase values and make continious (trend-like) data.
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

#constants------------------------------------------------------------------

timestep = 1.375


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


        intensity = np.append(intensity, np.sum(img))
           
intensity = remove_out_of_phase(intensity, 0.5)[0]

for i in range(len(time)):
        xticks.append(f"{time[i].minute}:{time[i].second}")

#Plot-figures-----------------------------------------------------------------

size = 16

plt.figure(figsize = (8, 4))
plt.title(f"Intensity profile 2021-05-22T06:15:08", size = size * 1.1, weight = "bold")
plt.plot(time, intensity, label = "full_sun")
plt.xticks(time, xticks)
plt.xlabel("time, [min:sec]", size = size)
plt.tick_params(axis = 'x', rotation = 45)
plt.ylabel("Intensity [arb. val.]", size = size)
plt.legend(fontsize = size * 0.5)
plt.tight_layout()
# plt.savefig(f"./SRH_I_{date}{frame}.png", transparent = False, dpi = 400, bbox_inches = "tight")
plt.show()