"""
Runs fullSWOF using a template located in folder `case_name`,
and saves output to a folder `out_dir`

Domain: 2D, with gaussian-filter generated blobs 
- Vary roughness contrast

Add in stripes again. 
Illustrate breakthrough.
Add code to calibrate to Cd 

"""
import os
import sys
import shutil
import itertools as it
import contextlib
import json
from multiprocessing import Pool

my_modules = ["write_SWOF"]
for mod in my_modules:
    if mod in sys.modules: 
        del sys.modules[mod]

from write_SWOF import *
from resistance_functions import *

# project directory
project_dir = "/Users/octaviacrompton/Dropbox/FullCSWOF/Tests"
out_name = "GRL_vary_2"
out_dir = os.path.join(project_dir, out_name)
 
default = {"L" : 50, 
           "l" : 200, 
           "dx" : 1,           
           "So" : 0.01,           
           "alpha_b" : "alpha_v",
           "alpha_v" : 24,
           "scale" : 1,                      
           "Ks_b" : 0.2,            
           "p" : 5,
           "sigma" : 3,
           "dt" : 60,
           "tr" : 30,
           "flux" : 1,
           "veg_type" : "blob",        
           "fV" : 0.2,        
           "fric" : 1,
           "t_rec" : 50,
           "D" : 0.005,           
           "scale" : 1,
           "scheme" : 'james',           
           "B_bc_init" : 2,           
           "aniso" : 1,
           "imax_scale" : 2,
           "phi_veg" : 0.1
          }

sim_dict = {"tr" : [10, 30, 60],
            "p" :  [1, 3, 5],             
            "Ks_v" : [3, 6], 
            "Ks_b" : [0.1, 0.5, 1], 
            "fV" : [0.1, 0.3, 0.5],
            "sigma" : [2, 3, 4],
            "aniso" : [1, 3],        
            "phi_veg" : [0.05, 0.15],
            "So" : [0.01]              
            }            

sim_vars = sim_dict.keys()

sim_list = [dict(list(zip(sim_vars, prod))) for prod in
                  it.product(*(sim_dict[var_name] 
                               for var_name in sim_vars))]


import random

# sim_list = random.sample(sim_list, 1000)
# print(len(sim_list))

# N = 500
# sim_list = {}
# for key, item in sim_dict.items():
#     low, high = item
#     sim_dict[key] = np.random.rand(100)*(high - low) + low
    
for i in range(len(sim_list)):
    sim_list[i]['seed'] = i

def call_SWOF(sdict):
    """
    wrapper function executes fullSWOF
    """      
    path = os.path.join(out_dir, "sim_list.json")
	
    with open(path, 'w') as fout:
        json.dump(sim_list , fout)

    path = os.path.join(out_dir, "default.json")
    with open(path, 'w') as fout:
        json.dump(default , fout)

    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])

    sim_dir = os.path.join(out_dir, sim_name)

    shutil.rmtree(sim_dir, ignore_errors=True)  
        
    os.mkdir (sim_dir) if not os.path.isdir(sim_dir) else 0
    if not os.path.isdir((sim_dir+ "/Outputs")):
        os.mkdir (sim_dir + "/Outputs") 
    if not os.path.isdir((sim_dir+ "/Inputs")):
        os.mkdir (sim_dir + "/Inputs") 
    
    # update parameter dictionary for this simulation case
    default.update(sdict) 
    
    if default['p'] < default['Ks_b']:
        print ("skip", sim_name)
        return

    if default['Ks_v'] < default['Ks_b']:
        print ("skip", sim_name)
        return        

    D = default["D"]
    phi_veg =  default["phi_veg"]
    So =  default["So"]
    
    if 'p-Ks_v' in default:
        default['Ks_v'] = default['p'] - default['p-Ks_v']
        
    if default["scheme"] == "james":

        default["alpha_o"] = get_alpha_o(default["D"], default["phi_veg"])
        default["R_v"] = np.pi/4*(1-default["phi_veg"])*default["D"]/default["phi_veg"]
        default['alpha_v'] = 24  
        default["fric"] = 6       

    if default["alpha_b"] ==  "alpha_v":
        # print ("veg and bare soil areas converge to same limit")
        default["alpha_b"] =  default["alpha_v"]
    
    # default['seed'] = default['sigma'] 

    default['T'] = (default['tr']+ default['t_rec'])*60    
    p0 = Params(".", default, overwrite = 1, T = default['T']) 
    
    # update parameter dictionary for this simulation case
    x,y,z = write_planar_topo(p0, sim_dir)

    if p0.veg_type == "uniform":
        veg = write_veg_uniform(p0, sim_dir)
        

    elif p0.veg_type == "blob":
        veg = write_veg_blob_aniso(p0, sim_dir)

    # write Ksat file, Ks.txt
    Ks = write_Ks(p0, sim_dir, veg)

    # write roughness parameter file, alpha.txt
    alpha = write_alpha(p0, sim_dir, veg)

    # storm with constant instensity and fixed duration
    write_rain(p0, sim_dir)
    write_huv(p0, sim_dir)
    #write_BC_B(p0, sim_dir)
    #write_BC_line(p0, sim_dir)
    write_BC_closed(p0, sim_dir)    
    write_BC_stop(p0, sim_dir)

    if "dt" in p0.__dict__:
        default["nbtimes"] =  p0.T/p0.dt + 1  
    else:
        default["dt"] =  p0.T/(p0.nbtimes + 1)    

    default["Nxcell"] = p0.Nxcell
    default["Nycell"] = p0.Nycell       
    default["T"] = p0.T      
    default["imaxcoef"] = p0.Ks_v*p0.imax_scale/3.6e5 # maximum infiltration rate
    write_input(default, case_dir = ".", sim_dir = sim_dir)        

    start = time.time()
    exec_name = "/Users/octaviacrompton/Dropbox/FullCSWOF/bin/FullSWOF_2D" # ../../bin/FullSWOF_2D    
    cmdout = os.popen("cd {0}; {1} ".format(sim_dir, exec_name)).read()
    end = time.time()
    os.remove(os.path.join(sim_dir, "Outputs", "flux_boundaries_LR.dat"))
    os.remove(os.path.join(sim_dir, "Outputs", "flux_boundaries_BT.dat"))

    cmdout =  process_cmdout(cmdout), (end-start)
    print("finishing: " + sim_name)
    return  cmdout

def process_cmdout(cmdout):
  
    return [l for l in cmdout.split("\n") if 
            ("%] done" not in l) and (len(l)>0)]


if __name__ == "__main__":

    # txt = input("Delete previous bank? y/n ")
    txt = "y"
    if txt != "y":
        print ("Adding to existing files")
    else:
        # remove all prev simulations
        print ("Removing existing files")
        shutil.rmtree(out_dir, ignore_errors=True) 
        os.mkdir (out_dir) if not os.path.isdir(out_dir) else 0
    
    shutil.copy(__file__, out_dir)

    os.mkdir(out_dir + "/code")
    shutil.copy("write_SWOF.py", out_dir  + "/code")
    shutil.copy("read_SWOF.py", out_dir  + "/code")


    print("starting simulations")
    # call_SWOF(sim_list[0])
    start = time.time()
    from multiprocessing import Pool
    pool = Pool()
    res = pool.map(call_SWOF, sim_list)
    pool.close()
    dt = time.time() - start                  
    print (dt)