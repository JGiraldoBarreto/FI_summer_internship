"""Provide functions related to Cryo-BIFE"""

import sys
import numpy as np
from tqdm import tqdm
import scipy.optimize as so
import matplotlib.pyplot as plt
from typing import Callable, Tuple


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
        log_likelihood = np.sum(np.log(np.dot(prob_mat.T, rho)))
        log_prior = kappa * prior_fxn(fe_prof)

        neg_log_posterior = -(log_likelihood + log_prior)

        return neg_log_posterior

    @staticmethod
    def grad2(
            path: np.ndarray,
            Sample_per_nodes: np.ndarray,
            images: np.ndarray,
            fe_prof: np.ndarray,
            sigma: float,
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
        prob_matrix = CryoBife.likelihood(path, images, Sample_per_nodes, sigma)
        Norm = np.dot(np.exp(-G),np.mean(prob_matrix,axis=1))
        
        grad = np.zeros((number_of_nodes,2))
        
        for k in range(number_of_nodes):

            Exppected_1_1 = (1/prob_matrix.shape[1]) * np.dot(Cv_differences.T[0,:,k],prob_matrix[k])
            Exppected_1_2 = (1/prob_matrix.shape[1]) * np.dot(Cv_differences.T[1,:,k],prob_matrix[k])

            Exppected_2_1 = np.mean(prob_matrix,axis=1)[k,:]*np.mean(Cv_differences[k],axis=0)[0]
            Exppected_2_2 = np.mean(prob_matrix,axis=1)[k,:]*np.mean(Cv_differences[k],axis=0)[1]  

            gradx = 2 * nu * np.exp(-G[k]) * (Exppected_1_1 - Exppected_2_1)/Norm
            grady = 2 * nu * np.exp(-G[k]) * (Exppected_1_2 - Exppected_2_2)/Norm

            grad[k,0] = np.sum(gradx) 
            grad[k,1] = np.sum(grady) 
            
        return -grad


    def optimizer(
            self,
            path: np.ndarray,
            images: np.ndarray,
            Sample_per_nodes: np.ndarray,
            sigma: float,
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

        prob_mat = self.likelihood(path, images, Sample_per_nodes, sigma)
        
        prob_mat = np.mean(prob_mat,axis=1)

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

kappa_2 = float(sys.argv[1])

images = np.loadtxt(sys.argv[2])    
ini_path = np.loadtxt(sys.argv[3])

Inv_T = float(sys.argv[4])

first_path = np.copy(ini_path)
sim_path = np.copy(ini_path)

nu = 20
normal_deviation = 1/nu

sigma = 0.5
step_size = 0.001
sampling_steps = 151
gradient_steps = 3
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

for s in tqdm(range(gradient_steps)):

    Cv_differences = []
    Sample_per_nodes = []

    for j in range(number_of_nodes):

        cv_diff = []
        cv_sampling = []

        ini_cv=np.array([sim_path[j,0],sim_path[j,1]])
        X = np.random.normal(ini_cv[0],normal_deviation,sampling_steps)
        Y = np.random.normal(ini_cv[1],normal_deviation,sampling_steps)

        for i in range(sampling_steps):

            new_cv = np.array([X[i],Y[i]])
            diff = new_cv - ini_cv
            cv_sampling.append(new_cv)
            cv_diff.append(diff)

        Sample_per_nodes.append(cv_sampling)
        Cv_differences.append(cv_diff)

    Sample_per_nodes = np.array(Sample_per_nodes)
    Cv_differences = np.array(Cv_differences)

    fe_prof, log_posterior = CryoBife.optimizer(CB,sim_path, images, Sample_per_nodes, sigma)
    Free_Energy.append(fe_prof)
    Log_Post.append(-log_posterior)

    #BATCH

    for i in range(number_of_batches):

        images_batch = images_shuffled[i*batch_size:(i+1)*batch_size]
        Gra = CryoBife.grad2(sim_path, Sample_per_nodes, images_batch, fe_prof, sigma, nu, Cv_differences)
        e_distance, e_grad = dist_energy_and_grad(sim_path,kappa_2,0)
        Gra = Gra + e_grad
        sim_path += -gradient_step_size * Gra * mask

    if residual_batches != 0:

        images_batch = images_shuffled[(number_of_batches-1)*batch_size:]
        Gra = CryoBife.grad2(sim_path, Sample_per_nodes, images_batch, fe_prof, sigma, nu, Cv_differences)
        e_distance, e_grad = dist_energy_and_grad(sim_path,kappa_2,0)
        Gra = Gra + e_grad        
        sim_path += -gradient_step_size * Gra * mask

    #END BACTH

    #Gra = CryoBife.grad2(sim_path, Sample_per_nodes, images, fe_prof, sigma, nu, Cv_differences)
    #e_distance, e_grad = dist_energy_and_grad(sim_path,kappa_2,0)
    #Gra = Gra + e_grad    
    #sim_path += -gradient_step_size*Gra*mask
    SP = np.copy(sim_path)
    Paths.append(SP)

Paths = np.array(Paths)

fe_prof, log_posterior = CryoBife.optimizer(CB,SP, images, Sample_per_nodes, sigma)
Free_Energy.append(fe_prof)
Log_Post.append(-log_posterior)

Log_Post = np.array(Log_Post)
Free_Energy = np.array(Free_Energy)

np.save('Path_BIFE_CV-Gauss_sampling_kappa_%s_T_%s.npy' %(kappa_2,Inv_T), Paths)
np.savetxt('FE_Path_BIFE_CV-Gauss_sampling_kappa_%s_T_%s' %(kappa_2,Inv_T), Free_Energy)
np.savetxt('LogPost_Path_BIFE_CV-Gauss_sampling_kappa_%s_T_%s' %(kappa_2,Inv_T), Log_Post)
np.save('rmsd_md_path.txt', Paths[-1])

for n in range(number_of_nodes):

    f = open('New_restraint_files/node_%s_rmsd' %(n+1),'w')
    f.write(str(Paths[-1,n,0])+' '+str(Paths[-1,n,1]))
    f.close()
