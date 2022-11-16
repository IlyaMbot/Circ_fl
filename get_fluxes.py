import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator as interpolate

images = np.load("cutted_images.npy")

image = np.copy(images[0])

x = np.array()
y = np.array()


'''
i = 5
fr_ch = i * 20

img1 = images[i]
img2 = images[i - 20]

img1[img1 > 50] = 50
img1[img1 < 40] = 0
img2[img2 > 50] = 50
img2[img2 < 40] = 0
img2 = np.roll(img2, shifts[i][0], axis = [0])
img2 = np.roll(img2, shifts[i][1], axis = [1])
print(np.sum(np.abs(img2-img1)))

plt.figure()
plt.imshow(img2 - img1)
plt.show()
'''