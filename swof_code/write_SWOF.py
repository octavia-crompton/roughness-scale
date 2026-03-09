from scipy.ndimage import gaussian_filter
import time
import os
import numpy as np
import re
import pandas as pd
import contextlib
np.random.seed(2)
# import cv2
from skimage.measure import block_reduce

from scipy.stats.mstats import gmean
flatten = lambda l: [item for sublist in l for item in sublist]

# from resistance_functions import *

class Params(object):
    """
    Class holds parameters for a SWOF simulation.
    
    Initialized using a parameter file located in a specified directory 
    (sim_dir).

    """
    def __init__(self, sim_dir, updates, overwrite = False, 
            mode = "writing", T = None):
        """
        Initialize with sim_dir containing parameter file
        """

        try:
            in_file = os.path.join(sim_dir , 'parameters.txt')
        except:
            in_file = os.path.join(sim_dir , 'parameters.dat')

        lines = read_lines(in_file)
        
        for i, line in enumerate(lines):
            if "::" in line:
                var_desc, var_name = line.split("::")[0].split("<")
                var_name = var_name.replace(">", "").strip()
                var_desc = var_desc.strip()

                val = line.split("::")[-1].strip()
                    
                if not val:
                    pass
                else:
                    try: 
                        val = float(val)
                    except ValueError:
                        pass

                setattr(self, var_name, val)
              
        for key  in updates:
            if overwrite:
                setattr(self, key, updates[key])    
            else:
                if key not in self.__dict__:
                    setattr(self, key, updates[key])    
        
        attrs = self.__dict__.keys()
        
        if "seed" not in self.__dict__:
            setattr(self, "seed", 0)                            

        if "sigma" not in self.__dict__:
            setattr(self, "sigma", 0)                                        
        
        if "dx" not in attrs:
            setattr(self, "dx", self.L/self.Nxcell)    
            
        if "dt" not in attrs:
            setattr(self, "dt", self.T/self.nbtimes)    

        self.Nxcell = int(np.round(self.L/self.dx))
        self.Nycell = int(np.round(self.l/self.dx)  )
        
        if mode == "writing":
            #self.T = (self.tr+self.t_rec)*60
            self.T = T

        sim_name = os.path.split(sim_dir)[0]
        sim_name = os.path.split(sim_name)[1]
        
        self.key = sim_name
        self.name = sim_name
        
    # def print(self):
        
    #     for attr, value in self.__dict__.items():
    #         if attr not in [ "doc", "topo"]:
    #             print(attr, '=', value)
             
    def update_dims(self, L, l, dx):

        self.L = L
        self.l = l        
        
        self.dx = dx
        
        self.Nxcell = int(self.L/self.dx)
        self.Nycell = int(self.l/self.dx)


def write_alpha(p0, sim_dir, veg):
    """
    Writes surface roughness file, alpha.txt

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory

    Notes: p0 contains 
        hillslope grid dimensions (Nxcell = width, Nycell = length)
        roughness values (alpha_b and alpha_v)
        pattern parameters (sigma = lengthscale, fV = veg fraction)
    """    
    x, y = make_grid(p0)    

    alpha = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 
    
    alpha[veg == 0] = p0.alpha_b
    alpha[veg == 1] = p0.alpha_v    

    alpha_array = np.vstack((x.ravel(), y.ravel(), alpha.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "alpha.txt")

    with open(filename, 'w') as f:
        f.write("# x y alpha[i][l] \n")
        
        for item in alpha_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
    
    # assert( np.abs(alpha.mean() - p0.alpha ) < 0.01)
    return alpha 


def write_veg_uniform(p0, sim_dir, zeros = False):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    

    np.random.seed(p0.seed)

    veg = np.ones((p0.Nxcell, p0.Nycell))
    
    if p0.fV == 0:
        # print ("fV = 0")
        veg = np.ones((p0.Nxcell, p0.Nycell))*0

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 

def write_veg(p0, sim_dir, veg):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    
    np.random.seed(p0.seed)

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 


def write_veg_blob(p0, sim_dir):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    

    np.random.seed(p0.seed)

    veg = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 
    gauss = gaussian_filter(veg, sigma=p0.sigma)
    threshold = np.percentile(gauss, p0.fV*100)
    veg[gauss > threshold] = 0
    veg[gauss <= threshold] = 1

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 

def write_veg_blob_aniso(p0, sim_dir):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    x, y = make_grid(p0)    

    np.random.seed(p0.seed)

    veg = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 
    r = abs(p0.aniso) if p0.aniso != 0 else 1
    if p0.aniso > 0:
        sig = (p0.sigma, p0.sigma / r)      # contour-aligned
    elif p0.aniso < 0:
        sig = (p0.sigma / r, p0.sigma)      # gradient-aligned
    else:
        sig = (p0.sigma, p0.sigma)           # isotropic
    gauss = gaussian_filter(veg, sigma=sig)
    threshold = np.percentile(gauss, p0.fV*100)
    veg[gauss > threshold] = 0
    veg[gauss <= threshold] = 1

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 


def write_veg_blob_buffer(p0, sim_dir):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    
    
    np.random.seed(p0.seed)

    veg = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 
    gauss = gaussian_filter(veg, sigma=p0.sigma)

    threshold = np.percentile(gauss, p0.fV*100)
    # wgt = 1 - (4*(p0.Nxcell ) + 4*(p0.Nycell))/(p0.Nxcell*p0.Nycell)
    # threshold = threshold/wgt

    veg[gauss >= threshold] = 0
    veg[gauss < threshold] = 1
    veg[:2, :] = 0
    veg[-2:, :] = 0
    veg[:, :2] = 0
    veg[:, -2:] = 0
    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 



def write_veg_dot(p0, sim_dir):
    """
    generates veg field and writes to input vegetation file, veg.txt
    

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)  
    np.random.seed(p0.seed)      

    from skimage import draw
    arr = np.zeros_like(x)
    rr, cc = draw.disk( 
            (int(x.shape[0]/2), int(y.shape[1]/2) + p0.offset), radius=int(p0.r/p0.dx), shape=arr.shape)
    arr[rr, cc] = 1
    veg = arr

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 

def write_veg_rectangle(p0, sim_dir):
    """
    generates veg field and writes to input vegetation file, veg.txt
    
    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)  
    np.random.seed(p0.seed)      

    from skimage import draw
    arr = np.zeros_like(x)

    w = int(p0.width/p0.dx)
    l = int(p0.length/p0.dx)
    # rr, cc = draw.rectangle((int(x.shape[0]/2 - w/2) ,  int(x.shape[1]/2 - l) ), 
    #                extent=(w, l), shape=x.shape)    
    rr, cc = draw.rectangle((int(x.shape[0]/2 - w/2) ,  int(x.shape[1]/2 -1 ) ), 
                   extent=(w, l), shape=x.shape)

    arr[rr, cc] = 1
    veg = arr

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 

def write_crescent(p0, sim_dir):
    """
    generates footprint field and writes to input footprintetation file, footprint.txt
    
    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)  
    np.random.seed(p0.seed)      

    from skimage import draw
    import cv2
    arr = np.zeros_like(x)

    width = p0.width
    rr, cc = draw.disk( 
            (int(x.shape[0]/2), int(y.shape[1]/2) - int(p0.radius/p0.dx) ), 
            radius=int(p0.radius/p0.dx), 
                shape=arr.shape)
    
    arr[rr, cc] = 1

    kernel = np.array([[0, 1, 0],
                      [1, 1, 1],
                      [0, 1, 0]], dtype = np.uint8)

    footprint = arr - cv2.erode(arr, np.ones((3,3)))
    footprint[:, :y.shape[1]//2 - int(width/p0.dx/2) ] = 0
    footprint[:int(x.shape[0]/2) - int(width/p0.dx/2), : ] = 0
    footprint[int(x.shape[0]/2) + int(width/p0.dx/2)+1:, : ] = 0 

    footprint_array = np.vstack((x.ravel(), y.ravel(), footprint.ravel())).T

    filename = os.path.join(sim_dir, "Inputs", "footprint.txt")

    with open(filename, 'w') as f:
        f.write("# x y berm[i][l] \n")
        
        for item in footprint_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     

    return footprint 


def write_veg_checkerboard(p0, sim_dir):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    x, y = make_grid(p0)  
    np.random.seed(p0.seed)      

    veg = checkerboard_scale(p0.Nxcell, p0.Nycell, p0.nchecks, p0.patch_scale )
    
    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T

    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg 



def checkerboard(Nxcell, Nycell, nchecks = 8, patch_length = 4):
    """

    """
    # patch_length =  int(np.sqrt(p0.fV*p0.Nxcell*p0.Nycell/p0.nchecks**2))
    ssize = int(Nxcell/nchecks)
    template = np.zeros([ssize, ssize])
    center = int(np.floor(ssize/2))
    template[center, center] = 1

    kernel = np.ones([patch_length,patch_length], dtype =np.uint8)

    tilex = int(Nxcell/ssize/2)
    tiley =  int(Nycell/ssize/2)
    template = cv2.dilate(template, kernel)
    template = np.hstack((template, 0*template))
    template = np.vstack((template, np.fliplr(template)))
    template = (np.tile(template, [tilex, tiley]))
    
    return template

def checkerboard_scale(Nxcell, Nycell, nchecks = 8, patch_scale = 0.5):
    """
    
    """
    ssize = int(Nxcell/nchecks)
    patch_length = int(ssize*patch_scale/2)
    template = np.zeros([ssize, ssize])
    center = int(np.floor(ssize/2))
    template[center, center] = 1
    template[center - patch_length: center+ patch_length+ 1, 
                center - patch_length:center+patch_length+1] = 1

    tilex = int(Nxcell/ssize/2)
    tiley =  int(Nycell/ssize/2)
    #template = cv2.dilate(template, kernel)
    template = np.hstack((template, 0*template))
    template = np.vstack((template, np.fliplr(template)))
    template = (np.tile(template, [tilex, tiley]))
    
    return template


def write_veg_randv(p0, sim_dir):
    """
    generates veg field and writes to input file, Ks.txt. 
    
    veg is generated using a Gaussian filter, with lengthscale sigma.

    Similar to write_veg_blob, but matches randv in my SVE model

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    
    np.random.seed(p0.seed)
    veg = sp.rand(p0.Nxcell+1, p0.Nycell+1) >= 1 - p0.fV
    veg = veg.astype(float)    

    gauss = gaussian_filter(veg.astype(float),
                                  sigma=(p0.sigma, p0.sigma))

    veg = (gauss > np.percentile(gauss, 100 * (1 - p0.fV))).astype(int)[:-1, :-1]


    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
     
    return veg.astype(float)


def write_veg_h_stripe(p0, sim_dir):
    """
    Writes surface roughness file, veg.txt
    """
    x, y = make_grid(p0)    

    # initialize veg sampling from uniform distribution
    veg = np.zeros(y.shape)

    stripe_count = int(p0.stripe_count)
    downslope = p0.downslope
    
    Nycell = p0.Nycell
    fV = p0.fV 
    
    veg_width = int(Nycell * fV / stripe_count  )

    bare_width = int(Nycell * (1. - fV) / stripe_count)
    distance_between_stripes = Nycell / stripe_count
    
    try:
        offset = p0.offset
    except:
        offset = 0
    
    if downslope == 'bare':
        lower_limit = bare_width

    elif downslope == 'veg':
        lower_limit = 0
    
    else:
        print ("specify a valid vegetation orientation")
        return

    start = Nycell -1 - lower_limit + offset
    
    if start < 0:
        start = 0
            
    for ind in np.arange(veg_width):

        inds = np.arange(start - ind, 0, # lower_limit + offset, 
            - distance_between_stripes, dtype=int)
        
        veg[:, inds] = 1

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
    
    return veg.astype(float)


def write_veg_v_stripe(p0, sim_dir):
    """
    Writes surface roughness file, veg.txt
    """
    x, y = make_grid(p0)    

    # initialize veg sampling from uniform distribution
    veg = np.zeros(y.shape)

    stripe_count = int(p0.stripe_count)
    downslope = p0.downslope
    
    Nxcell = p0.Nxcell
    fV = p0.fV 
    
    veg_width = int(Nxcell * fV / stripe_count)
    bare_width = int(Nxcell * (1. - fV) / stripe_count)
    distance_between_stripes = Nxcell / stripe_count

    if downslope == 'bare':
        lower_limit = 0
    elif downslope == 'veg':
        lower_limit = bare_width
    else:
        print ("specify a valid vegetation orientation")
        return

    for ind in np.arange(veg_width):
        inds = np.arange(lower_limit + ind, Nxcell, 
            distance_between_stripes, dtype=int)
        veg[inds, :] = 1


    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
    
    return veg  


def write_Ks(p0, sim_dir, veg):
    """
    Writes Ks to input file, Ks.txt, using input veg map
    
    Ks is generated using a Gaussian filter, with lengthscale sigma.

    Inputs:
        p0: Params object
        sim_dir : SWOF simulation directory
    """
    # 
    x, y = make_grid(p0)    
    
    Ks = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 

    Ks[veg < 1] = p0.Ks_b/3.6e5
    Ks[veg == 1] = p0.Ks_v/3.6e5    
    Ks_array = np.vstack((x.ravel(), y.ravel(), Ks.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "Ks.txt")

    with open(filename, 'w') as f:
        f.write("# x y Ks[i][l] \n")
        
        for item in Ks_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
    
    # assert( np.abs(Ks.mean() - p0.Ks ) < 0.01)
    return Ks 

def write_planar_topo(p0, sim_dir):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)
    z = rescale(1-y)*p0.So*p0.l
    

    topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    filename = os.path.join(sim_dir, "Inputs", "topography.txt")

    with open(filename, 'w') as f:
        f.write("# x y z[i][l] \n")
        for item in topo_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x,y,z

def write_wall_topo(p0, sim_dir, veg):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)
    z = rescale(1-y)*p0.So*p0.l
    z[veg == 1] = z[veg == 1] + p0.wall_height

    import cv2
    smooth = p0.smooth
    z = cv2.GaussianBlur(z,(smooth,smooth),0)

    center = int(p0.Nxcell/2)
    if p0.gap == 1:
        z[center] = rescale(1-y[center])*p0.So*p0.l
    
    topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    filename = os.path.join(sim_dir, "Inputs", "topography.txt")

    with open(filename, 'w') as f:
        f.write("# x y z[i][l] \n")
        for item in topo_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x,y,z    

def make_curvilinear_berm(p0, footprint, wall_height = 2.5):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)

    L=p0.l # length of the hill slope
    xmax=p0.L
    slope = p0.So
    H = slope*L
    a = p0.a

    omega = - a*H/(2*L) # curvature paramter 
        
    # z = H*(1-y/L)+omega*x**2
    y_edge = y # - p0.dx/2
    x_edge = x # - p0.dx/2
    z = H*(1-y_edge/L)+omega*( x_edge - xmax/2)**2
        
    ####  adding a berm to Dana's topography
    import cv2

    if wall_height > 0:
        # minimum elevation in berm footprint plus berm height
        try:            
            minh = z[footprint == 1].min() + wall_height
            z[footprint == 1] = minh
        except:
            print ("fail", sim_dir.split("/")[-1])
            return
    
    smooth = p0.smooth
    z = cv2.GaussianBlur(z,(smooth,smooth),1)
    z = z.round(5)

    return x,y,z


def write_curvilinear_berm(p0, sim_dir, footprint):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)
    # x = x - p0.dx/2
    # y = y - p0.dx/2

    L=p0.l # length of the hill slope
    xmax=p0.L
    slope = p0.So
    H = slope*L
    a = p0.a

    omega = - a*H/(2*L) # curvature paramter 
        
    # z = H*(1-y/L)+omega*x**2
    y_edge = y # - p0.dx/2
    x_edge = x # - p0.dx/2
    z = H*(1-y_edge/L)+omega*( x_edge - xmax/2)**2
        
    ####  adding a berm to Dana's topography
    import cv2

    if p0.wall_height > 0:
        # minimum elevation in berm footprint plus berm height
        try:            
            minh = z[footprint == 1].min() + p0.wall_height
            z[footprint == 1] = minh
        except:
            print ("fail", sim_dir.split("/")[-1])
            return

    smooth = p0.smooth
    z = cv2.GaussianBlur(z,(smooth,smooth),1)
    z = z.round(5)

    topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    filename = os.path.join(sim_dir, "Inputs", "topography.txt")

    with open(filename, 'w') as f:
        f.write("# x y z[i][l] \n")
        for item in topo_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x,y,z

def write_curvilinear_topo(p0, sim_dir = ''):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)
    # x = x - p0.dx/2
    # y = y - p0.dx/2

    L= p0.l # length of the hill slope
    xmax=p0.L
    slope = p0.So
    H = slope*L
    a = p0.a

    omega= - a*H/(2*L) # curvature paramter 
    y_edge = y # - p0.dx/2
    x_edge = x # - p0.dx/2
    z = H*(1-y_edge/L)+omega*( x_edge - xmax/2)**2

    if sim_dir != '':
        topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

        filename = os.path.join(sim_dir, "Inputs", "topography.txt")

        with open(filename, 'w') as f:
            f.write("# x y z[i][l] \n")
            for item in topo_array:
                f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x,y,z
 

def write_topo(p0, sim_dir, gauss):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    from read_SWOF import get_patchL
    x, y = make_grid(p0)
    z = rescale(1-y)*p0.So*p0.l
    z = z + gauss
    z = z.round(5)

    # z = fix_adversity(z)
    topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    filename = os.path.join(sim_dir, "Inputs", "topography.txt")

    with open(filename, 'w') as f:
        f.write("# x y z[i][l] \n")
        for item in topo_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x, y, z


def fix_adversity(z):
    """
    Potentially useful for microtopography with adverse slope conditions
    """
    bad = np.where((np.diff(z, axis = 1) > 0)[:, -1])[0]

    grad = np.diff(z, axis = 1) > 0

    patch_LB, upslope_V = get_patchL(grad, 1000)   
    adverse = patch_LB[:, -1]

    for ind in bad:
        z[ind, - int(adverse[ind])-1:]  =  z[ind, - int(adverse[ind])]
        # print (z[ind, - int(adverse[ind]):])

    grad = np.diff(z, axis = 1) > 0

    return  z


def write_shelf_topo(p0, sim_dir):
    """
    Write fullSWOF input coordinates (x,y,z) for a planar hillslope,
    for a simulation instance located in sim_dir
    """    
    x, y = make_grid(p0)
    z = rescale(1-y)*p0.So*p0.l
    z = z.round(5)
    z[:, :3] = z[:, 2:3]
    
    topo_array = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    filename = os.path.join(sim_dir, "Inputs", "topography.txt")

    with open(filename, 'w') as f:
        f.write("# x y z[i][l] \n")
        for item in topo_array:
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))
    
    return x,y,z

    
def rescale (zed):
    """ 
    rescales an input array from 0 to 1
    """ 
    return (zed- zed.min())/(zed.max() - zed.min())

def write_huv(p0, sim_dir):
    """
    Writes input h,u,v file, huv.txt, 
    with all variables initialized to 0 
    """
    x, y = make_grid(p0)    
    if "h0" in p0.__dict__:
        h0 = p0.h0
    else:
        h0 = 0
            
    h = np.ones_like(y)*h0
    u = np.ones_like(y)*0
    v = np.ones_like(y)*0

    huv_array = np.vstack((x.ravel(), y.ravel(), h.ravel(), u.ravel(), v.ravel())).T

    filename = os.path.join(sim_dir, "Inputs", "huv.txt")
    with open(filename, 'w') as f:
        f.write("#x y h[i][l] u[i][l] v[i][l]\n")
        for item in huv_array:
            f.write("%s %s %s %s %s \n" % (item[0], item[1], item[2], item[3], item[4]))

def write_huv_dam(p0, sim_dir):
    """
    Writes input h,u,v file, huv.txt, 
    with all variables initialized to 0 
    """
    x, y = make_grid(p0)    

    h = np.ones_like(y)*0
    u = np.ones_like(y)*0
    v = np.ones_like(y)*0
    h[int(p0.Nxcell/2)-1:int(p0.Nxcell/2)+2, 1] =  p0.h
    huv_array = np.vstack((x.ravel(), y.ravel(), h.ravel(), u.ravel(), v.ravel())).T

    filename = os.path.join(sim_dir, "Inputs", "huv.txt")
    with open(filename, 'w') as f:
        f.write("#x y h[i][l] u[i][l] v[i][l]\n")
        for item in huv_array:
            f.write("%s %s %s %s %s \n" % (item[0], item[1], item[2], item[3], item[4]))

def make_grid_copy(dx, Nxcell, Nycell):
    """
    creates a grid of cell centers
    """     
    x = np.linspace(dx/2, (Nxcell-0.5)*dx, Nxcell) 
    y = np.linspace(dx/2, (Nycell-0.5)*dx, Nycell) 
    
    # x = np.linspace(dx/2, p.L - dx, Nxcell) 
    # y = np.linspace(dx/2, p.l - dx, Nycell)     

    x, y = np.meshgrid(x, y)
    x = x.T
    y = y.T

    x = np.round(x, 3)
    y = np.round(y, 3)

    return x, y

def make_grid(p):
    """
    creates a grid of cell centers
    """     
    x = np.linspace(p.dx/2, (p.Nxcell-0.5)*p.dx, int(p.Nxcell) )
    y = np.linspace(p.dx/2, (p.Nycell-0.5)*p.dx, int(p.Nycell) )
    
    # x = np.linspace(p.dx/2, p.L - p.dx, p.Nxcell) 
    # y = np.linspace(p.dx/2, p.l - p.dx, p.Nycell)     

    x, y = np.meshgrid(x, y)
    x = x.T
    y = y.T

    x = np.round(x, 3)
    y = np.round(y, 3)

    return x, y


def write_single_stripe(p0, sim_dir):
    """
    Writes surface roughness file, veg.txt
    """
    x, y = make_grid(p0)    
    
    # initialize veg sampling from uniform distribution
    veg = np.zeros(y.shape)
    
    Nycell = p0.Nycell
    fV = p0.fV 
    
    veg_width = int(Nycell * fV )

    try:
        offset = p0.offset
    
    except:
        offset = 0
    
    start = Nycell - 1 -  int(offset/p0.dx)
    
    if start < 0:
        start = 0
            
    inds = np.arange(start , start - veg_width , - 1, dtype=int)
        
    veg[:, inds] = 1

    veg_array = np.vstack((x.ravel(), y.ravel(), veg.ravel())).T
    
    filename = os.path.join(sim_dir, "Inputs", "veg.txt")

    with open(filename, 'w') as f:
        
        f.write("# x y veg[i][l] \n")
        
        for item in veg_array:
            
            f.write("%s %s %s \n" % (item[0], item[1], item[2]))    
    
    return veg.astype(float)

def write_rain(p0, sim_dir):
    """
    writes rain.txt input file for SWOF,
    rainfall has constant intensity p_mps (m/s) for duration T (s)
    """
    tr_s = p0.tr*60
    p_mps = p0.p/3.6e5       
    T = p0.T
    
    filename = os.path.join(sim_dir, "Inputs", "rain.txt")
    with open(filename, 'w') as f:
            f.write("%s %s \n" % (0, p_mps))    
            f.write("%s %s \n" % (tr_s - 1, p_mps))                
            f.write("%s %s \n" % (tr_s, 0 ))                            
            f.write("# %s %s \n" % (T, 0 ))                                        


def write_rain_csv(rain_file, sim_dir):
    """
    writes rain.txt input file for SWOF,
    rainfall has constant intensity p_mps (m/s) for duration T (s)
    """    
    precip = pd.read_csv(os.path.join(rain_file))[['Time', 'Precipitation']]  
    

    T = precip.Time.max() + 10
    filename = os.path.join(sim_dir,  "Inputs", "rain.txt")
    with open(filename, 'w') as f:
        for index, row in precip.iterrows():

            t, p = row

            f.write("%s %s \n" % (t*60, p/1000/3600))                                          


        f.write("# %s %s \n" % (T*60, 0 ))                                                                                   
    return precip

def write_input(sdict, case_dir, sim_dir):   
    """
    modify template parameter file located in case_dir, and save to sim_dir

    sdict : dictionary of distinct parameters
    case_dir : directory with template file (e.g., '.')
    sim_dir : directory of new simulation. 

    """
    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])
    
    sim_vars = sdict.keys()
    in_file =  case_dir + '/parameters.txt'
    out_file = sim_dir  + '/Inputs/parameters.txt'    
    
    # remove any existing parameter file
    with contextlib.suppress(FileNotFoundError):
        os.remove(out_file)  
    
    # loop over in_file lines and overwrite any files in out_file
    with open(in_file) as f:
        with open(out_file, "w") as f1:
            for line in f:
                for sim_var in sim_vars:
                    test_string =  "<"+ sim_var + ">"
                    if test_string in line:
                                            
                        split0  = line.split("::")[0]
                        line = ":: ".join([split0, str(sdict[sim_var])]) + "\n"

                f1.write(line)

                
def read_lines(infile):
    """
    Returns infile as a list of lines 
    """
    fin = open(infile) 
    lines = fin.readlines()
    fin.close()
    return lines


############################
# INHOMOGENEOUS BOUNDARIES

def write_BC_B(p0, sim_dir):
    """
    writes BC_B.txt input file for SWOF,
    """
    t_inflow = p0.t_inflow

    filename = os.path.join(sim_dir, "Inputs", "BC_B.txt")
    with open(filename, 'w') as f:
            f.write("0.0      BC_B_1.txt \n")    
            f.write("{0}      BC_B_2.txt \n".format(t_inflow*60))                


def write_BC_closed(p0, sim_dir):
    """
    writes BC_B.txt input file for SWOF, closed
    """

    filename = os.path.join(sim_dir, "Inputs", "BC_B.txt")
    with open(filename, 'w') as f:
            f.write("0.0      BC_B_2.txt \n")                

     

def write_BC_start(p0, sim_dir  ):
    """
    """
    xc, yc = make_grid(p0)  
    xc_bound = xc[:, 1]
    BC = np.ones_like(xc_bound, dtype = int)*2
    q = np.zeros_like(xc_bound, dtype = int)
    width = int(p0.width/2)
    
    mid = int(len(xc_bound)/2)
    BC[mid - width : mid + width] = 5
    bottom_imp_discharge = p0.bottom_imp_discharge
    
    h_normal = p0.bottom_imp_h
    filename = os.path.join(sim_dir, "Inputs", "BC_B_1.txt".format(1))
    
    with open(filename, 'w') as f:
            f.write("#x      c      q   h \n")    
            for i, xc in enumerate(xc_bound):
                if BC[i] == 2:
                    f.write("{0}     {1}       \n".format(xc, BC[i]))    
                else:
                    f.write("{0}     {1}      {2}      {3} \n".format(xc, BC[i], bottom_imp_discharge, h_normal)  )  
                     

def write_BC_line(p0, sim_dir):
    """
    bottom_imp_discharge = inflow
    Units in m3/s (per grid cell)

    To get m3/s (per boundary) – inflow/dx*p0.L
    """
    xc, yc = make_grid(p0)  
    xc_bound = xc[:, 1]
    BC = np.ones_like(xc_bound, dtype = int)
    q = np.zeros_like(xc_bound, dtype = int)
    
    BC[:] = 5
    bottom_imp_discharge = p0.bottom_imp_discharge
    
    h_normal = p0.bottom_imp_h
    filename = os.path.join(sim_dir, "Inputs", "BC_B_1.txt".format(1))
    
    with open(filename, 'w') as f:
            f.write("#x      c      q   h \n")    
            for i, xc in enumerate(xc_bound):
                f.write("{0}     {1}      {2}      {3} \n".format(xc, BC[i], bottom_imp_discharge, h_normal)  )  
                     
             
        
def write_BC_stop(p0,  sim_dir , i = 2 ):
    """
    """
    xc, yc = make_grid(p0)  
    xc_bound = xc[:, 1]
 
    filename = os.path.join(sim_dir, "Inputs", "BC_B_{0}.txt".format(i))
    
    with open(filename, 'w') as f:
            f.write("#x      c      q   h \n")    
            for i, xc in enumerate(xc_bound):
                
                f.write("{0}     2        \n".format(xc))    
                
def write_TC_B(p0, sim_dir):
    """
    writes TC_B_T.txt input file for SWOF – 
        helps take care of adverse BCs at the outlet
        needed for microtopography
    """
    t_inflow = p0.t_inflow

    filename = os.path.join(sim_dir, "Inputs", "TC_B.txt")
    with open(filename, 'w') as f:
            f.write("0.0      TC_B_1.txt \n")                
               
def write_TC_start(p0, sim_dir, topo):
    """
    writes TC_B_1.txt input file for SWOF – 
        helps take care of adverse BCs at the outlet
    """
    from read_SWOF import get_patchL
    xc, yc = make_grid(p0)  
    xc_bound = xc[:, -1]
    
    grad = (topo[:, 1:] - topo[:,  :-1]) > 0

    patch_LB, upslope_V = get_patchL(grad, 1000)   
    adverse = patch_LB[:, -1]

    filename = os.path.join(sim_dir, "Inputs", "TC_B_1.txt".format(1))
    
    with open(filename, 'w') as f:
            f.write("#x      c      q   h \n")    
            for i, xc in enumerate(xc_bound):
                if adverse[i] > 0:
                    f.write("{0}     {1}     \n".format(xc, 2))    
                else:
                    f.write("{0}     {1}     \n".format(xc, 3 )) 
                     
    
