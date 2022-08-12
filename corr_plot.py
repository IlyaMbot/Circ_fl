from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os, glob
import antlib

filename = glob.glob('./fluxes/2021_05_22*.fits')
filename = sorted(filename, key=os.path.basename)

with fits.open(filename[0], memmap = False) as f:
    f.verify('silentfix')
    data = f[0].data

time = data[0]
data = np.sum(np.abs(data), 0)

plt.figure()
plt.plot(time, data)
plt.show()