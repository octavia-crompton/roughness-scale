'''
Runs fullSWOF using a template located in folder `case_name`,
and saves output to a folder `out_dir`

Domain: 2D, with gaussian-filter generated blobs 
- Vary roughness contrast

Add in stripes again. 
Illustrate breakthrough.
Add code to calibrate to Cd 

'''
import os
import sys
import shutil
import itertools as it
import contextlib
import json
from multiprocessing import Pool

my_modules = ['write_SWOF']
for mod in my_modules:
    if mod in sys.modules: 
        del sys.modules[mod]

from write_SWOF import *
from resistance_functions import *

# project directory
project_dir = '/Users/octaviacrompton/Dropbox/FullCSWOF/Tests'
out_name = 'micro_So-2'
out_dir = os.path.join(project_dir, out_name)
 
default = {'L' : 100, 
           'l' : 100, 
           'dx' : 2,           
           'So' : 0.01,           
           'alpha_v' : 0.05,
           'scale' : 1,                      
           'Ks_b' : 0.,  
           'Ks_v' : 5,                     
           'dt' : 30,
           'flux' : 1,
           'rec' : 1,
           'veg_type' : 'dot',        
           't_inflow' : 0, # may want to add a storm... 
           'fric' : 1,
           'D' : 0.005,           
           'scale' : 1,           
           'Ks_v' : 3,
           't_rec' : 50,                      
           'aniso' : 1,
           'sigma' : 2, 
           'tr' : 30,
           'p' : 1,
           'fV' : 0. ,
           'seed' : 0,           
           'T_bc_init' : 2, # constant
           'B_bc_init' : 2  # constant 
          }


# vegetation fraction 
sim_dict = {
            "alpha_b" : [0.03, 0.05, 0.07],
            "radius" : [5, 10],
            "p" : [2, 4, 6, 8],   
            "Ks_b" : [0.5, 1], 
            "So" : [0.02, 0.005],
            "sigma" : [2, 5],
            "amplitude" : [0., 0.1, 0.2]        
            }


sim_vars = sim_dict.keys()

sim_list = [dict(list(zip(sim_vars, prod))) for prod in
                  it.product(*(sim_dict[var_name] 
                               for var_name in sim_vars))]

# for i in range(len(sim_list)):
#     sim_list[i]['seed'] = i

def call_SWOF(sdict):
    '''
    wrapper function executes fullSWOF

    '''       

    path = os.path.join(out_dir, 'sim_list.json')
	
    with open(path, 'w') as fout:
        json.dump(sim_list , fout)

    path = os.path.join(out_dir, 'default.json')
    with open(path, 'w') as fout:
        json.dump(default , fout)

    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])

    sim_dir = os.path.join(out_dir, sim_name)

    shutil.rmtree(sim_dir, ignore_errors=True)  
        
    os.mkdir (sim_dir) if not os.path.isdir(sim_dir) else 0
    if not os.path.isdir((sim_dir+ '/Outputs')):
        os.mkdir (sim_dir + '/Outputs') 
    if not os.path.isdir((sim_dir+ '/Inputs')):
        os.mkdir (sim_dir + '/Inputs') 
    
    # update parameter dictionary for this simulation case
    case = default.copy()
    case.update(sdict) 

    D = case['D']    
    So =  case['So']
    
    if 'p-Ks_v' in case:
        case['Ks_v'] = case['p'] - case['p-Ks_v']

    case['T'] = (case['tr']+ case['t_rec'])*60

    if case['So'] <  1e-4:
        case['Lbound'] = 3
        case['Rbound'] = 3
        case['Bbound'] = 3

    p0 = Params(".", case, overwrite = 1,  T = case['T'])  

    # update parameter dictionary for this simulation case

    
    if p0.veg_type == 'blob':

        veg = np.random.uniform(0, 1, (p0.Nxcell, p0.Nycell)) 

        bound = 5
        veg[:bound] = veg.mean()
        veg[-bound:] = veg.mean()
        veg[:, :bound] = veg.mean()
        veg[:, -bound:]= veg.mean()

        gauss = gaussian_filter(veg, sigma=(p0.sigma, p0.sigma/p0.aniso))
        threshold = np.percentile(gauss, (1-p0.fV)*100)
        
        topo = (gauss - gauss.mean())*p0.amplitude/gauss.std()
        veg[gauss > threshold] = 1
        veg[gauss <= threshold] = 0

        x,y,z = write_topo(p0, sim_dir, topo)

        veg = write_veg(p0, sim_dir, veg)


    elif p0.veg_type == 'dot':

        x, y = make_grid(p0)    

        arr = np.zeros_like(x)
        x_center, y_center = int(p0.dx*x.shape[0]/2),  int(p0.dx*y.shape[1]/2)

        a = p0.radius*p0.dx
        c = p0.amplitude

        z = (1 - ((x-x_center)**2)/a**2 - ((y-y_center)**2)/a**2)
        z[z < 0] = 0
        z = np.sqrt(z)
        
        veg = z > 0

        z = c*(1 - ((x-x_center)**2)/a**2 - ((y-y_center)**2)/a**2)
        z[z < 0] = 0
        z = np.sqrt(z)

        z = gaussian_filter(z, sigma=p0.sigma)

        # z = c*(z - z.min())/(z.max() - z.min())
        topo = z
        
        x,y,z = write_topo(p0, sim_dir, topo)

        veg = write_veg(p0, sim_dir, veg)
        print (veg.mean(), topo.max(), topo.min(),
             np.round(np.diff(topo, axis = 1).max(), 2))


    # write Ksat file, Ks.txt
    Ks = write_Ks(p0, sim_dir, veg)

    # write roughness parameter file, alpha.txt
    alpha = write_alpha(p0, sim_dir, veg)

    # storm with constant instensity and fixed duration
    write_rain(p0, sim_dir)
    write_huv(p0, sim_dir)
    write_BC_B(p0, sim_dir)
    write_BC_line(p0, sim_dir)
    write_BC_stop(p0, sim_dir)
    write_TC_B(p0, sim_dir)
    write_TC_start(p0, sim_dir, topo)
    
    if 'dt' in p0.__dict__:
        case['nbtimes'] =  p0.T/p0.dt + 1  
    else:
        case['dt'] =  p0.T/(p0.nbtimes + 1)    

    case['Nxcell'] = p0.Nxcell
    case['Nycell'] = p0.Nycell       
    case['T'] = p0.T      
    case['imaxcoef'] = p0.Ks_v/3.6e5 # maximum infiltration rate
    write_input(case, case_dir = '.', sim_dir = sim_dir)        

    start = time.time()
    exec_name = '/Users/octaviacrompton/Dropbox/FullCSWOF/bin/FullSWOF_2D' # ../../bin/FullSWOF_2D    
    cmdout = os.popen('cd {0}; {1} '.format(sim_dir, exec_name)).read()

    end = time.time()
    os.remove(os.path.join(sim_dir, 'Outputs', 'flux_boundaries_LR.dat'))
    os.remove(os.path.join(sim_dir, 'Outputs', 'flux_boundaries_BT.dat'))

    cmdout =  process_cmdout(cmdout), (end - start)
    print('finishing: ' + sim_name)
    return  cmdout


def process_cmdout(cmdout):
  
    return [l for l in cmdout.split('\n') if 
            ('%] done' not in l) and (len(l)>0)]


if __name__ == '__main__':

    # txt = input('Delete previous bank? y/n ')
    txt = 'y'
    if txt != 'y':
        print ('Adding to existing files')
    else:
        # remove all prev simulations
        print ('Removing existing files')
        shutil.rmtree(out_dir, ignore_errors=True) 
        os.mkdir (out_dir) if not os.path.isdir(out_dir) else 0
    
    shutil.copy(__file__, out_dir)

    os.mkdir(out_dir + '/code')
    os.mkdir(out_dir + '/animations')
    os.mkdir(out_dir + '/static')

    shutil.copy('write_SWOF.py', out_dir  + '/code')
    shutil.copy('read_SWOF.py', out_dir  + '/code')

    print('starting simulations')
    # call_SWOF(sim_list[0])
    start = time.time()
    from multiprocessing import Pool
    pool = Pool()
    res = pool.map(call_SWOF, sim_list)
    pool.close()
    dt = time.time() - start                  
    print (dt)