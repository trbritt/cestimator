from scipy.io import matlab
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

def two_dim_ellipsoid(location, square_dispersion, scale, axis, color, label, plot_eig_vectors=False, plot_square=False):
   """
   Calculates and plots a 2D ellipsoid with optional features.

   Args:
       location (np.ndarray): Center location of the ellipsoid.
       square_dispersion (np.ndarray): Square root of the dispersion matrix.
       scale (float): Scaling factor for the ellipsoid.
       plot_eig_vectors (bool, optional): Flag to plot eigenvectors (default: False).
       plot_square (bool, optional): Flag to plot a square around the ellipsoid (default: False).
   """

   # Calculate eigenvalues and eigenvectors
   eigenvalues, eigenvectors = np.linalg.eigh(square_dispersion)
   # Define angles for parametric representation
   angles = np.linspace(0, 2*np.pi, 501)

   # Create centered ellipse points
   centered_ellipse = np.zeros((len(angles), 2))
   for i in range(len(angles)):
       y = np.array([np.cos(angles[i]), np.sin(angles[i])])
       centered_ellipse[i] = eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ y

   # print(centered_ellipse)
   # Calculate actual ellipse points with location and scaling
   ellipse_points = location + scale * centered_ellipse

   # Plot the ellipse
   axis.plot(ellipse_points[:, 0], ellipse_points[:, 1], color=color, linewidth=2, label=label)
   # Plot the square (if desired)
   if plot_square:
       dispersion = np.sqrt(np.diag(square_dispersion))
       vertex_lr = location + scale * dispersion
       vertex_ll = location - scale * dispersion
       vertex_ul = location + np.array([scale * dispersion[0], -scale * dispersion[1]])
       vertex_ur = location + np.array([-scale * dispersion[0], scale * dispersion[1]])
       square = np.vstack([vertex_lr, vertex_ll, vertex_ul, vertex_ur, vertex_lr])
       axis.plot(square[:, 0], square[:, 1], color=color, linewidth=2)

   # # Plot eigenvectors (if desired)
   if plot_eig_vectors:
       l1 = scale * np.sqrt(eigenvalues[0])
       l2 = scale * np.sqrt(eigenvalues[1])
       sign = np.sign(eigenvectors[0, 0])  # Handle potential reflection

       # Eigenvector 1
       start_a, end_a = location[0], location[0] + sign * eigenvectors[0, 0] * l1
       start_b, end_b = location[1], location[1] + sign * eigenvectors[0, 1] * l1
       axis.plot([start_a, end_a], [start_b, end_b], color=color, linewidth=2)

       # Eigenvector 2
       start_a, end_a = location[0], location[0] + eigenvectors[2, 0] * l2
       start_b, end_b = location[1], location[1] + eigenvectors[2, 1] * l2
       axis.plot([start_a, end_a], [start_b, end_b], color=color, linewidth=2)

   # Set aspect ratio and labels
   axis.set_aspect('equal')
   axis.set_xlabel('X')
   axis.set_ylabel('Y')

   # # Add legend (if applicable)
   if plot_square or plot_eig_vectors:
       axis.legend()

def hubertM(x):
   """
   Calculates the location  and scatter matrix for multivariate distribution using Huber-Hampel outlier penalizer
   
   Args:
       X (np.ndarray): Input data matrix (T rows, N columns).

   Returns:
       tuple: A tuple containing the shrunk mean vector (M_Shr) and shrunk covariance matrix (S_Shr).
       
   """
   def OutlierCutoff(d, d0):
       index = np.where(d<=d0, True, False)
       b_2 = 1.25
       omega = (d0/d) * np.exp(-0.5*((d-d0)**2/(b_2**2)))
       omega[index] = 1
       return omega
   tolerance = 10e-6
   error = 1e6
   N, T = X.shape
   w = np.ones((T,))
   Zeros = np.zeros((N,))
   mu = Zeros 
   sigma = np.zeros((N,N))
   
   d0 = np.sqrt(N) + np.sqrt(2)
   
   while error>tolerance:
       
       mu_old = mu 
       sigma_old = sigma 
       
       mu = Zeros
       for t in range(T):
           mu += w[t] * x[:,t].T
       mu /= w.sum()
       
       sigma = np.zeros((N,N)) 
       for t in range(T):
           tmp = (x[:,t].T - mu).reshape(2,1)
           sigma += w[t] * w[t] * tmp @ tmp.T
       sigma /= (w.T @ w)
       inv_sigma = np.linalg.inv(sigma)
       d = []
       for t in range(T):
           d.append(np.sqrt((x[:,t].T-mu).T @ inv_sigma @ (x[:,t].T-mu)))
       w = OutlierCutoff(d, d0)
       error = np.trace(np.linalg.matrix_power(sigma-sigma_old, 2) + (mu-mu_old).T @ (mu-mu_old))
   return mu, sigma

def shrinkage(x):
   """
   Calculates shrinkage estimates for mean and covariance matrix.

   Args:
       X (np.ndarray): Input data matrix (T rows, N columns).

   Returns:
       tuple: A tuple containing the shrunk mean vector (M_Shr) and shrunk covariance matrix (S_Shr).
   """

   N, T = X.shape  # Get dimensions of the data
   # Calculate mean without centering
   mu_no_par = np.mean(x, axis=1)

   # Calculate unregularized covariance matrix
   std_no_par = np.cov(x)

   # Target vector (all zeros)
   b = np.zeros((N,))

   # Eigenvalues of the unregularized covariance
   lambda_hat = np.linalg.eigvals(std_no_par)
   # # Calculate weight for mean shrinkage
   a = (1/T) * (lambda_hat.sum() - 2*lambda_hat.max()) / ( mu_no_par.T @ mu_no_par )
   a = max(0,min(a,1))
   mu_shr = (1-a)*mu_no_par + a * b
   C = np.mean(lambda_hat)*np.eye(N)
   #compute optimal weight
   num = 0
   for t in range(T):
       num += 1/T* np.trace(np.linalg.matrix_power(x[:,t] @ x[:,t].T - std_no_par, 2))
   denom = np.trace((std_no_par - C) @ (std_no_par - C))
   a = (1/T) * num / denom
   a = max(0, min(a,1))
   std_shr = (1-a)*std_no_par + a*C
   # C = n
   return mu_shr, std_shr

def MLE(x):
   """
   this function computes the maximum likelihood estimate of the dof Nu, location mu
   and scatter matrix sigma of T iid observations from the n-variate t distribution
   """
   def recursion(x, nu, tolerance):
       N, T = x.shape
       w = np.ones((T,1))
       mu = np.zeros((N,1))
       sigma = np.zeros((N,N))
       error = 1e6
       while error > tolerance:
           mu_old = mu
           sigma_old = sigma 
           
           W = w @np.ones((1,N))
           mu = (np.sum(W * x.T, axis=0).T / np.sum(w)).reshape(2,1)
           x_c = x - np.transpose(np.ones((T,1))@mu.T )
           sigma = (W * x_c.T).T @ x_c.T / T 
           
           inv_sigma = np.linalg.inv(sigma)
           ma2 = np.sum((x_c.T@inv_sigma)*x_c.T, axis=1)
           w = (nu + N) / (nu + ma2)
           error = np.trace(np.linalg.matrix_power(sigma-sigma_old, 2)/N + (mu-mu_old).T @ (mu-mu_old)/N)
       return mu, sigma 
   
   def log_likelihood(x, nu, mu, sigma):
       inv_sigma = np.linalg.inv(sigma)
       N, T = x.shape
       norm = -N/2 * np.log(nu*np.pi) + np.log(gamma((nu+N)/2)) - np.log(gamma(nu/2)) - 0.5*np.log(np.linalg.det(sigma))
       retval = 0
       for t in range(T):
           centered = x[:,t].reshape(N,1) - mu
           ma2 = centered.T @ inv_sigma @ centered
           retval += norm - (nu+N)/2 * np.log(1+ma2/nu)
       return retval
   tolerance = 0.01*np.mean(np.percentile(x, 75)-np.percentile(x, 25))
   nus = [1,2,4,7,12,20]
   
   LL = []
   for idn, nu in enumerate(nus):
       m, s = recursion(x, nu, tolerance)
       LL.append(log_likelihood(x, nu, m, s))
   print(LL)
   index = np.argmax(LL)
   nu = nus[index]
   m, sigma = recursion(x, nu, tolerance)
   sigma *= nu / (nu-2)
   return m.T, sigma

data = np.loadtxt('../../examples/US_SwapRates.csv')
# lets choose which sets of rates we'll compare
choose_rates = [1, 3] #0+1=2yr, 1+1=5yr, 3+1=10yr
Y = np.vstack((data[:,choose_rates[0]], data[:,choose_rates[1]]))
X = Y[:,1:] - Y[:,:-1]

# First, we do non-parametric estimation of the relationship
mu_no_par = np.mean(X, axis=1)
std_no_par = np.cov(X)
# Second, we do shrinkage methods
mu_shr, std_shr = shrinkage(X)
mu_h, std_h = hubertM(X)
mu_MLE, std_MLE = MLE(X)

fig, ax = plt.subplots(1,1)
ax.scatter(X[0,:], X[1,:], marker='x', color='k', s=5)

two_dim_ellipsoid(mu_no_par, std_no_par, 2, ax, 'b', 'non-param', plot_eig_vectors=False, plot_square=False)
two_dim_ellipsoid(mu_shr, std_shr, 2, ax, 'g', 'shrinkage', plot_eig_vectors=False, plot_square=False)
two_dim_ellipsoid(mu_h, std_h, 2, ax, 'r', 'hubert', plot_eig_vectors=False, plot_square=False)
two_dim_ellipsoid(mu_MLE, std_MLE, 2, ax, 'y', 'MLE', plot_eig_vectors=False, plot_square=False)

ax.legend()
plt.show()