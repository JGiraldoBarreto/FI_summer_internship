# Run as python3 optimize_Gaussian_LogPost.py

import numpy as np
import numpy.random as r
from numpy.linalg import inv


def Gaussian_Images(Path, Num_Images_Matrix, factor, sigma):
    """
    Images with gaussian noise from the toy system model proposed by ALex.
    Parameters
    ----------
    Path : numpy array
           size L*2 numpy matrix of the nodes of a path with length L, in the coordinate space.
    Num_Images_Matrix : numpy array
           size N*M numpy matrix of the number of images per model of a system, with N*M the total number of models.
    sigma : float
           standard deviation of the gaussian model.
    factor : float
           Controls the final number of desired images (float>0)
    Returns:
    --------
    Gauss_Image_vector : numpy array
           size I*2 numpy matrix of the gaussian images I, in the coordinate space.
    """
    amount_of_images = 0
    Gauss_Image_vector = []
    for i in range(Path.shape[0]):
        Index0 = np.int(Path[i, 0])
        Index1 = np.int(Path[i, 1])
        amount_of_images = np.int(factor*Num_Images_Matrix[Index0][Index1])
        for j in range(amount_of_images):
            x_coord = Index0 + r.normal(0, sigma)
            y_coord = Index1 + r.normal(0, sigma)
            Gauss_Image_vector.append([x_coord, y_coord])
    Gauss_Image_vector = np.array(Gauss_Image_vector)
    return(Gauss_Image_vector)


def sample_grid_data(grid,Num_Images_Matrix,factor,sigma,N_total,inverse_T=None):

    # Generate_grid_data

    Grid = np.array(grid) #np.loadtxt("data/grid")  # Gridpoints of model
    Grid = Grid - np.ones(Grid.shape)  # correction to match python notation.
    Num_images = Num_Images_Matrix #np.loadtxt("data/Num_Images_grid-3well")

    if inverse_T is not None:
        Num_images = Num_images.astype('float') ** inverse_T
    Num_images *= (N_total / np.sum(Num_images))
    Num_images = np.round(Num_images)

    G_Vec = Gaussian_Images(Grid, Num_images, factor, sigma)

    return G_Vec, Num_images


if __name__ == "__main__":
    factor = 1.0  # scaling factor
    sigma = 1.0   # standard deviation
    Path = np.random.randint(20, size=(14, 2))  # dummy path
    Num_Images_Matrix = np.random.randint(50, size=(20, 20))  # dummy Matrix with number of images per model
    Gaussian_Images_Vector = Gaussian_Images(Path, Num_Images_Matrix, factor, sigma)
    print('Gaussian Images vector: \n', Gaussian_Images_Vector)
