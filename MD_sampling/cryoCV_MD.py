"""Provide functions related to Cryo-BIFE"""

import sys
import numpy as np
from tqdm import tqdm
from typing import Callable, Tuple

import MDAnalysis as mda
from MDAnalysis.analysis import align

import scipy.optimize as so
from scipy.spatial.transform import Rotation


class CryoBife:
    """CryoBife provides cryo-bife's prior, likelihood, posterior and
    the optimizer as described in 10.1038/s41598-021-92621-1. """

    @staticmethod
    def integrated_prior(fe_prof: np.ndarray) -> float:
        """Calculate the value of the prior for the given free-energy profile.

        :param fe_prof: Array with the values of the free-energy profile in each
            node of the path.

        :returns: The value of the prior for the free-energy profile.
        """
        acc_der_fe_prof = sum(np.diff(fe_prof)**2)
        log_prior = np.log(1 / acc_der_fe_prof**2)

        return log_prior

    @staticmethod
    def likelihood(
            path: np.ndarray,
            images: np.ndarray,
            Sample_per_nodes: np.ndarray, 
            sigma: float) -> np.ndarray:
        """Calculate cryo-bife's likelihood matrix given a path and a dataset of images

        :param path: Array with the values of the variables at each node of the path.
                     Shape must be (n_models, n_dimensions).
        :param images: Array with all the experimental images.
                       Shape must be (n_images, image_dimensions)
        :param sigma: Overall noise among the images

        :returns: Array with the likelihood of observing each image given each model.
                  Shape will be (n_images, n_models)
        """
        
        norm = 1 / (2 * np.pi * sigma**2)
        number_of_images = images.shape[0]
        number_of_nodes = path.shape[0]
        number_of_samples = Sample_per_nodes.shape[1]
        prob_matrix = np.zeros((number_of_nodes,number_of_samples,number_of_images))

        for i in range(number_of_nodes):
            for j in range(number_of_samples):

                prob_matrix[i,j] = norm * np.exp( (-0.5*1/(sigma**2))*\
                                                  ((Sample_per_nodes[i,j,0] - images[:,0])**2+\
                                                   (Sample_per_nodes[i,j,1] - images[:,1])**2))        

        return prob_matrix

    def neg_log_posterior(
            self,
            fe_prof: np.ndarray,
            kappa: float,
            prob_mat: np.ndarray,
            beta: float = 1,
            prior_fxn: Callable = None) -> float:
        """Calculate cryo-bife's negative log-posterior.

        :param fe_prof: Array with the values of the free-energy profile (FEP)
                        in each node of the path.
        :param beta: Temperature.
        :param kappa: Scaling factor for the prior.
        :prior_fxn: Function used to calculate the FEP's prior

        :returns: Value of the negative log-posterior
        """

        if prior_fxn is None:
            # Default is the prior from the paper
            prior_fxn = self.integrated_prior

        # TODO: Think for a better name for rho
        rho = np.exp(-beta * fe_prof) #density vec
        rho = rho / np.sum(rho) #normalize, Eq.(8)

        # Sum here since iid images; logsumexp
        log_likelihood = np.sum(np.log(np.dot(np.exp(prob_mat), rho)))
        log_prior = kappa * prior_fxn(fe_prof)

        neg_log_posterior = -(log_likelihood + log_prior)

        return neg_log_posterior

    @staticmethod
    def grad2(
            prob_matrix: np.ndarray,
            fe_prof: np.ndarray,
            nu: float,
            Cv_differences: np.ndarray,
            prior_fxn: Callable = None) -> np.ndarray:
        """Calculate gradient of cryo-bife's log-posterior following eq 16 of cryoBIFE_MD.

        :param path: Array with the values of the variables at each node of the path.
                     Shape must be (n_models, n_dimensions).
        :param fe_prof: Array with the values of the free-energy profile (FEP)
                        in each node of the path.
        :param images: Array with all the experimental images.
                       Shape must be (n_images, image_dimensions)
        :param sigma: standar deviationof the synthetic images.
        :prior_fxn: Function used to calculate the FEP's prior

        :returns: Array with the values of the gradient.
        """
        G = fe_prof
        number_of_nodes = G.shape[0]
        Norm = np.dot(np.exp(-G),np.mean(prob_matrix,axis=0))
        
        grad = np.zeros((number_of_nodes,2))
        
        for k in range(number_of_nodes):

            Exppected_1_1 = (1/prob_matrix.shape[0]) * np.dot(Cv_differences[0,:,k],prob_matrix[:,k,:])
            Exppected_1_2 = (1/prob_matrix.shape[0]) * np.dot(Cv_differences[1,:,k],prob_matrix[:,k,:])

            Exppected_2_1 = np.mean(prob_matrix,axis=0)[k,:]*np.mean(Cv_differences[:,:,k],axis=1)[0]
            Exppected_2_2 = np.mean(prob_matrix,axis=0)[k,:]*np.mean(Cv_differences[:,:,k],axis=1)[1]  

            gradx = 2 * nu * np.exp(-G[k]) * (Exppected_1_1 - Exppected_2_1)/Norm
            grady = 2 * nu * np.exp(-G[k]) * (Exppected_1_2 - Exppected_2_2)/Norm

            grad[k,0] = np.sum(gradx) 
            grad[k,1] = np.sum(grady) 
            
        return -grad


    def optimizer(
            self,
            path: np.ndarray,
            prob_mat: np.ndarray,
            initial_fe_prof: np.ndarray = None) -> np.ndarray:
        """Find the optimal free-energy profile given a path and a dataset of images

        :param path: Array with the values of the variables at each node of the path.
                     Shape must be (n_models, n_dimensions).
        :param images: Array with all the experimental images.
                       Shape must be (n_images, image_dimensions)
        :param sigma: TODO.
        :param fe_prof: Initial guess for the free-energy profile

        :returns: Optimized free-energy profile
        """

        kappa = 1
        n_models = path.shape[0]

        if initial_fe_prof is None:

            initial_fe_prof = 1.0 * np.random.randn(n_models)
        
        prob_mat = np.mean(prob_mat,axis=2)

        optimized_fe_prof = so.minimize(self.neg_log_posterior,
                                        initial_fe_prof,
                                        method='CG',
                                        args=(kappa, prob_mat))

        return (optimized_fe_prof.x,optimized_fe_prof.fun)


def dist_energy_and_grad(
        path: np.ndarray,
        kappa_2: float,
        equib_dist: float) -> Tuple[float, np.ndarray]:
    """Calculate harmonic distance constraint energy and grad for a path.

    :param path: Array with the initial values of the free-energy profile in each
                node of the path.
    :kappa_2: harmonic constant for the constraint
    :equib_dist: equilibrium distance between nodes of the path

    :returns: harmonic distance constraint energy and gradient
    """

    # Calculate energy
    energy_dist = np.sum((np.sqrt(np.sum((path[:-1] - path[1:])**2, axis=1)) - equib_dist)**2)
    energy_dist *= 0.5 * kappa_2

    # Calculate distances between two or more nodes
    def distance(x1, x2):
        return np.sqrt(np.sum((x1 - x2)**2, axis=-1))

    # Calculate gradient
    grad_dist = np.zeros_like(path)

    grad_dist[0] = (1 - equib_dist / distance(path[0], path[1])) *\
                (path[0] - path[1])

    grad_dist[-1] = (1 - equib_dist / distance(path[-1], path[-2])) *\
                    (path[-1] - path[-2])

    grad_dist[1:-1] = (1 - equib_dist / distance(path[1:-1], path[2:])[:,None]) *\
                    (path[1:-1] - path[2:]) +\
                    (1 - equib_dist / distance(path[1:-1], path[:-2])[:,None]) *\
                    (path[1:-1] - path[:-2])

    grad_dist *= kappa_2

    return energy_dist, grad_dist


########### MAIN ############

#### Part 1: Setting the MD data###

#### 1.1: Aligning models

ref_model = mda.Universe("total_traj_frames/traj_1.pdb")
ref_model.atoms.translate(-ref_model.atoms.center_of_mass())

n_pdbs = 5689
models = np.empty((n_pdbs, 3, len(ref_model.atoms)))

for i in range(0, n_pdbs):

    tmp_model = mda.Universe(f"total_traj_frames/traj_{i+1}.pdb")
    align.alignto(tmp_model, ref_model, select="all", weights="mass")

    models[i] = tmp_model.select_atoms("all").positions.T

#### 1.2 Defining functions to generate images


def gen_quat(num_quaternions):
    #Sonya's code
    
    quaternions = np.zeros((num_quaternions, 4))
    count = 0

    while count < num_quaternions:

        quat = np.random.uniform(-1, 1, 4) #note this is a half-open interval, so 1 is not included but -1 is
        norm = np.sqrt(np.sum(quat**2))

        if ( 0.2 <= norm <= 1.0 ):
            quaternions[count] = quat/norm
            count += 1

    return quaternions

def calc_ctf(n_pixels, pixel_size, amp, phase, b_factor):

    freq_pix_1d = np.fft.fftfreq(n_pixels, d=pixel_size)
    freq_x, freq_y = np.meshgrid(freq_pix_1d, freq_pix_1d)

    freq2_2d = freq_x**2 + freq_y**2

    env = np.exp(-b_factor * freq2_2d * 0.5)

    ctf = amp * np.cos(phase * freq2_2d * 0.5) -\
          np.sqrt(1 - amp**2) * np.sin(phase * freq2_2d * 0.5) + np.zeros_like(freq2_2d) * 1j
    
    return ctf * env #/ amp

def apply_ctf(image, ctf):

      image_ft = np.fft.fft2(image)
      image_ft *= ctf

      image_ctf = np.fft.ifft2(image_ft).real

      return image_ctf

def add_noise(image, snr):

    mean_image = np.mean(image)
    std_image = np.std(image)

    mask = np.logical_or(
        image >= mean_image + 0.5 * std_image, image <= mean_image - 0.5 * std_image
    )

    signal_mean = np.mean(image[mask])
    signal_std = np.std(image[mask])

    noise_std = signal_std / np.sqrt(snr)
    noise = np.random.normal(loc=signal_mean, scale=noise_std, size=image.shape)

    img_noise = image + noise

    return img_noise

def create_image(coord, quat, n_pixels, pixel_size, sigma):

    coord_copy = coord.copy()

    rot_matrix = np.array(Rotation.from_quat(quat).as_matrix())

    coord_copy = np.matmul(rot_matrix, coord_copy)

    norm = 1 / (2 * np.pi * sigma**2 * coord_copy.shape[1])

    grid_min = -pixel_size * (n_pixels - 1) * 0.5
    grid_max = pixel_size * (n_pixels - 1) * 0.5 + pixel_size

    grid = np.arange(grid_min, grid_max, pixel_size)

    gauss_x = np.exp(-0.5 * (grid[:, None] - coord_copy[0, :]) ** 2 / sigma**2)
    gauss_y = np.exp(-0.5 * (grid[:, None] - coord_copy[1, :]) ** 2 / sigma**2)

    image = np.matmul(gauss_x, gauss_y.T) * norm

    return image

#### 1.3 Defining imaging parameters and checking that everything works

# Image generator
N_PIXELS = 128
PIXEL_SIZE = 0.25
SIGMA_IMGS = 0.7

# Noise
SNR = 0.1

# CTF
AMP = 0.1
B_FACTOR = 0.0
DEFOCUS = 0.1
ELECWAVEL = 0.019866
PHASE = DEFOCUS * np.pi * 2.0 * 10000 * ELECWAVEL

quat = gen_quat(1)[0]
image = create_image(models[0], quat, N_PIXELS, PIXEL_SIZE, SIGMA_IMGS)

ctf = calc_ctf(N_PIXELS, PIXEL_SIZE, AMP, PHASE, B_FACTOR)
image_ctf = apply_ctf(image, ctf)

image_noise = add_noise(image, SNR)
# image_noise = add_noise(image_ctf, SNR) # CTF + NOISE

#### 1.4 Generating an image for every model with one random rotation

images = np.zeros((models.shape[0], N_PIXELS, N_PIXELS))
quats = gen_quat(models.shape[0])

ctf = calc_ctf(N_PIXELS, PIXEL_SIZE, AMP, PHASE, B_FACTOR)

for i in range(models.shape[0]):

    images[i] = create_image(models[i], quats[i], N_PIXELS, PIXEL_SIZE, SIGMA_IMGS)
    # images[i] = apply_ctf(images[i], ctf)
    images[i] = add_noise(images[i], SNR)

def calc_log_lklhood(
    path_samples, ref_images, ref_quats, n_pixels, pixel_size, sigma, cb_sigma
):

    norm = 1 / (2 * np.pi * cb_sigma**2)

    log_lkl = np.zeros(
        (ref_images.shape[0], path_samples.shape[0], path_samples.shape[1])
    )

    print('Generating Images')

    for i in tqdm(range(ref_images.shape[0])):
        for j in range(path_samples.shape[0]):
            for k in range(path_samples.shape[1]):

                tmp_image = create_image(
                    path_samples[j, k], ref_quats[i], n_pixels, pixel_size, sigma
                )

                log_lkl[i, j, k] = np.sum((tmp_image - ref_images[i]) ** 2) * (
                    -0.5 / cb_sigma**2
                ) + np.log(norm)

    return log_lkl

N_NODES = 10
N_SAMPLES = 200
CV_DIM = 2
CB_SIGMA = 0.5
SPRING_CONSTANT = 10 # nu

path_samples = np.empty((N_NODES, N_SAMPLES, 3, len(ref_model.atoms)))

for i in range(0, N_NODES):
    for j in range(0, N_SAMPLES):

        if i == 0:
            path_samples[i, j] = models[0]

        elif i == N_NODES - 1:
            path_samples[i, j] = models[-1]

        else:
            #tmp_model = mda.Universe(f"samples_per_nodes/node-{i+1}_samples/traj_frames/traj_{101 + j}.pdb")
            tmp_model = mda.Universe(f"node-{i+1}_folder/traj_frames/traj_{101 + j}.pdb")
            align.alignto(tmp_model, ref_model, select="all", weights="mass")

            path_samples[i,j] = tmp_model.select_atoms("all").positions.T

log_lklhood = calc_log_lklhood(path_samples, images[:15], quats[:15], N_PIXELS, PIXEL_SIZE, SIGMA_IMGS, CB_SIGMA).T

########################
########################
########################

kappa_2 = float(sys.argv[1])
ini_path = np.loadtxt(sys.argv[2])    

first_path = np.copy(ini_path)
sim_path = np.copy(ini_path)

nu = 20
normal_deviation = 1/nu

sigma = 0.5
step_size = 0.001
sampling_steps = 151
gradient_steps = 20
gradient_step_size = 0.0001
number_of_nodes = ini_path.shape[0]

#BATCH

batch_size = int(images.shape[0] * 0.1)
number_of_batches = images.shape[0]//batch_size
residual_batches = images.shape[0]%batch_size

images_shuffled = images.copy()
np.random.shuffle(images_shuffled)

#END BATCH

Paths = []
Log_Post = []
Free_Energy = []

mask = np.ones_like(ini_path)
mask[0] = np.zeros((2,))
mask[-1] = np.zeros((2,))

CB = CryoBife()

print('Optimization process:')


Cv_differences = np.zeros((CV_DIM,N_SAMPLES,N_NODES))

for n in range(1,number_of_nodes-1):

    #Cv_diff = np.loadtxt('node-%s_folder/Cv_differences_node-%s' %(n,n)).T #np.random.randn(CV_DIM,N_SAMPLES,N_NODES)
    #Cv_differences.apend(Cv_diff)
    Cv_differences[:,:,n] = np.loadtxt('../MD_sampling/node-%s_folder/Cv_differences_node-%s' %(n+1,n+1)).T[:,:N_SAMPLES]

Cv_differences = np.array(Cv_differences) ### ---> Estos son los valores de RMSD que da Plumed

for s in tqdm(range(gradient_steps)):

    fe_prof, log_posterior = CryoBife.optimizer(CB,sim_path, log_lklhood)
    Free_Energy.append(fe_prof)
    Log_Post.append(-log_posterior)

    Gra = CryoBife.grad2(log_lklhood, fe_prof, nu, Cv_differences)
    e_distance, e_grad = dist_energy_and_grad(sim_path,kappa_2,0)
    Gra = Gra + e_grad    
    sim_path += -gradient_step_size*Gra*mask
    SP = np.copy(sim_path)
    Paths.append(SP)

Paths = np.array(Paths)

fe_prof, log_posterior = CryoBife.optimizer(CB,SP, log_lklhood)
Free_Energy.append(fe_prof)
Log_Post.append(-log_posterior)

Log_Post = np.array(Log_Post)
Free_Energy = np.array(Free_Energy)

np.save('Path_BIFE_CV-Gauss_sampling_kappa_%s.npy' %(kappa_2), Paths)
np.savetxt('FE_Path_BIFE_CV-Gauss_sampling_kappa_%s' %(kappa_2), Free_Energy)
np.savetxt('LogPost_Path_BIFE_CV-Gauss_sampling_kappa_%s' %(kappa_2), Log_Post)
np.savetxt('rmsd_md_path.txt', Paths[-1])

for n in range(number_of_nodes):

    f = open('New_restraint_files/node_%s_rmsd' %(n+1),'w')
    f.write(str(Paths[-1,n,0])+' '+str(Paths[-1,n,1]))
    f.close()
