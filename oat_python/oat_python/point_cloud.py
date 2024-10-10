
import copy
from typing import Tuple
import numpy as np
import sklearn
from sklearn.neighbors import radius_neighbors_graph

import scipy




def sample_points_from_disk(radius, num_points):
    """
    Randomly samples points from the disk of given radius centered at the origin in the Euclidean plane.

    :param radius: the radius of the disk
    :param num_points: number of points to sample
    """
    # Generate random angles
    angles = np.random.uniform(0, 2*np.pi, num_points)
    
    # Generate random radii within the disk
    radii = np.sqrt(np.random.uniform(0, radius**2, num_points))
    
    # Convert polar coordinates to Cartesian coordinates
    x = radii * np.cos(angles)
    y = radii * np.sin(angles)
    
    return x, y


def two_circles():
    """
    Points evenly spaced along two tangent circles in the Euclidean plane.
    """
    theta  = np.linspace(0, 2*np.pi, 100)
    cloud = np.zeros((100,2))
    cloud[:,0] = np.cos(theta)
    cloud[:,1] = np.sin(theta)
    cloud2 = copy.deepcopy(cloud)
    cloud2[:,0] += 1
    cloud = np.concatenate((cloud,cloud2),axis=0)
    return cloud

def torus_curve( npoints=100 , radius_minor=0.3, radius_major=1, angle_initial=0, nturns=1, ):
    """
    Return a sequence of points that curves around a torus `nturns` times.
    """
    r           = radius_minor
    R           = radius_major
    T           =   np.linspace(0, 2 * np.pi,          npoints) # major angle
    t           =   np.linspace(0, 2 * np.pi * nturns, npoints) + angle_initial * np.pi # minor angle
    x           =   ( R + r * np.cos( t ) ) * np.cos( T )
    y           =   ( R + r * np.cos( t ) ) * np.sin( T )
    z           =   r * np.sin( t )
    cloud      =   np.zeros((npoints,3))
    cloud[:,0]=x; cloud[:,1]=y; cloud[:,2]=z;
    return cloud


def torus(radius_outer, radius_inner, npoints_outer, npoints_inner, repeat_last):
    """
    Returns numpy arrays x, y, z, each having size 
    - `(npoints_inner+1) x (npoints_inner+1)` (if `repeat_last = True`)
    - `(npoints_inner) x (npoints_inner)` (if `repeat_last = False`)
    For each `p`,
    - `(x[p], y[p], z[p])` is a circle that goes around the torus in the "long direction"
    - `(x[:,p], y[:,p], z[:,p])` is a circle that goes around the torus in the "short direction"    

    :param radius_outer: the larger radius
    :param radius_inner: the smaller radius
    :param npoints_outer: number of points to sample from the larger circle
    :param npoints_inner: number of points to sample from the smaller circle
    :param repeat_last: determines whether the first and last rows (respectively, columns) of each `x, y, z` matrix are equal; defualt = `True`
    """
    phi_values = np.linspace(0, 2 * np.pi, npoints_outer, endpoint=False)
    theta_values = np.linspace(0, 2 * np.pi, npoints_inner, endpoint=False)
    if repeat_last:
        phi_values = np.append(phi_values, 0)
        theta_values = np.append(theta_values, 0)
    phi, theta = np.meshgrid(phi_values, theta_values)

    x = (radius_outer + radius_inner * np.cos(theta)) * np.cos(phi)
    y = (radius_outer + radius_inner * np.cos(theta)) * np.sin(phi)
    z = radius_inner * np.sin(theta)

    return x, y, z

def spiral_sphere( npoints=50, embedding_dim=3, noise_scale=1, xaxisthresh=1, random_seed=None ):
    """
    A sample of points fro a sphere intersected with a halfspace of form `x ≤ xaxisthresh`.

    :param `npoints` int: number of points in the cloud
    :param `embedding_dim`: embedding dimension of the cloud, defaults to 3
    :param `noise_scale` float: points are generated in a deterministic fashion, then motified by addint a small amount of noise drawn in an iid fashion from the unit cube of size `[0,noise_scale]^n`.
    :param `xaxisthresh` float: determines how far on the x axis you go; xthresh=0 returns a hemisphere
    :param `random_seed`: a random seed for the generation of noise, if desired
    """

    if not random_seed is None:
        np.random.seed( random_seed )

    theta   = np.linspace(0, 9 * 2 * np.pi, npoints)
    x       = np.linspace(-1, xaxisthresh, npoints)
    cloud = np.zeros((npoints, embedding_dim))
    cloud[:,0] += x
    cloud[:,1] += np.cos(theta) * (1-x**2)
    cloud[:,2] += np.sin(theta) * (1-x**2)
    cloud += np.random.rand(npoints, embedding_dim) * noise_scale
    
    return cloud
    


def half_dome( npoints=50, noise_scale=1 ):
    """
    Samples points evenly from the half of the Euclidean 2-sphere that lies between x=-1 and x=0.

    This function is just a wrapper for `slice_of_sphere( npoints=npoints, xmin=-1, xmax=0) + noise.
    :param npoints: number of points to sample
    """
    cloud = slice_of_sphere( npoints=npoints, randomize=True, xmin=-1.0, xmax=0.0)
    cloud += np.random.rand( npoints, 3 ) * noise_scale
    return cloud   


# A function to generate random points on a sphere, sliced between the plane x=xmin and x=xmax
def slice_of_sphere(  npoints=1, xmin=-1.0, xmax=1.0, randomize=True, random_seed=None, ):
    """
    Samples points evenly from the portion of unit sphere in Euclidean 3-space that lies between the planes `x=xmin` and `x=xmax`, using the Fibonacci method.
    :param npoints: number of points to sample
    :param randomize: adds some randomization
    :xmin: (approximate)  lower threshold for x coordinates
    :xmax: (approximate)  upper threshold for x coordinates
    :randomize: randomizes an initial condition
    :random_seed: sets a random seed for initialization
    """
    rnd = 1.
    if randomize:
        if not random_seed is None:
            np.random.seed( random_seed )
        rnd = np.random.random() * npoints

    points = []
    increment_y   = (xmax - xmin)/npoints
    increment_phi = np.pi * (3. - np.sqrt(5.))

    for i in range( npoints ):
        x = xmin + i * increment_y # ((i * offset) - 1) + (offset / 2) #  #
        r = np.sqrt(1 - x**2)
        phi = ((i + rnd) % npoints) * increment_phi
        y = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append([x, y, z])

    return np.array(points)        



def annulus( npoints=1, rad0=1, rad1=2, random_seed=None ):
    """
    Sample points uniformly from an annulus
    """
    cloud = np.zeros((npoints,2))

    if not random_seed is None:
        np.random.seed( random_seed )

    scale_parameter     =   rad0**2-rad1**2
    for i in range(npoints):
        theta       =   2 * np.pi * np.random.rand()
        rnd         =   np.random.rand()
        rho         =   np.sqrt( rnd * scale_parameter + rad1**2 );
        cloud[i][0]    =   rho * np.cos(theta)
        cloud[i][1]    =   rho * np.sin(theta)

    return cloud

      