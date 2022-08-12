from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os, glob

#Read-I-and-V-files-----------------------------------------------------------

# filenames = glob.glob('./sandbox*/pictures/061415/*02_?.fits')
filenames = glob.glob('./sandbox*/pictures/061653/*07_?.fits')

filenames = sorted(filenames, key=os.path.basename)

with fits.open(filenames[0], memmap = True) as f:
        f.verify('silentfix')
        imgI = f[0].data

y, x = np.histogram(imgI, bins = 'auto')
x = x[:-1]

fond1 = np.exp(-x ** 2) 
fond2 = np.exp(-x ** 2 * 100) * 1000
fond3 = np.exp(-(x - 5 )** 2 * 50) * 2000






plt.figure()
plt.plot(x, y, "o")
plt.plot(x, fond1)
plt.plot(x, fond2)
plt.plot(x, fond3)

plt.show()