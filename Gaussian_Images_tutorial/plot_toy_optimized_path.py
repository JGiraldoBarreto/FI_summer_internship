import numpy as np
import matplotlib.pyplot as plt

Path = np.load('Path_BIFE_CV-Gauss_sampling_kappa_50.0_T_3.0.npy')
ini_path = np.loadtxt('data/O4')
Gaussian_Images = np.loadtxt('Gaussian_images_Inv_T_3.0')

plt.figure()

plt.hist2d(Gaussian_Images[:,0],Gaussian_Images[:,1],bins=(40,40))
plt.plot(ini_path[:,0],ini_path[:,1], 'k-*', label = 'Initial Path')
plt.plot(Path[-1,:,0],Path[-1,:,1], 'r-*', label = 'Optimized Path')

plt.colorbar()

plt.show()

