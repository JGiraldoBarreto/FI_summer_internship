import sys
import numpy as np
import numpy.random as r
from Generate_Gaussian_Images_Toymodel import sample_grid_data

sigma = 0.5
factor = 1.0
N_total=20000
Inv_Temp = float(sys.argv[3])
Full_grid = np.loadtxt(sys.argv[1])
Num_Images_Matrix = np.loadtxt(sys.argv[2])

Images, New_Image_Matrix = sample_grid_data(Full_grid,Num_Images_Matrix,factor,sigma,N_total,Inv_Temp) 

np.savetxt('Gaussian_images_Inv_T_%s' %Inv_Temp, Images)
np.savetxt('New_Matrix_Inv_T_%s' %Inv_Temp, New_Image_Matrix, fmt='%1.5f')

## If want to see the example distribution, please uncomment the following lines:

import matplotlib.pyplot as plt
Im = np.loadtxt('Gaussian_images_Inv_T_%s' %Inv_Temp)
plt.hist2d(Im[:,0],Im[:,1],bins=(40,40))
plt.colorbar()
plt.show()
