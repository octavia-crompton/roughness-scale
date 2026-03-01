"""
Runs fullSWOF using a template located in folder `case_name`,
and saves output to a folder `out_dir`

Domain: 2D, with gaussian-filter generated blobs 
- Vary roughness contrast

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
out_name = "scenarios_dot_manning_dx4"
out_dir = os.path.join(project_dir, out_name)
 
default = {"L" : 32.2, 
           "l" : 40.2, 
           "dx" : 0.4,           
           "So" : 0.01,           
           "alpha_b" : 0.03,         
           "Ks_b" : 0.5,            
           "Ks_v" : 5,
           "p" : 5,
           "offset" : 5,           
           "dt" : 15,
           "flux" : 1,
           "veg_type" : "dot",
           "seed" : 1,                    
           "fV" : 0.05,        
           "D" : 0.004, 
           "t_inflow" : 0,                     
           "phi_veg" :  0.1,
           "r" : 5,
           "fric" : 1,
           'B_bc_init' : 2,
           "t_rec" :  15                  
          }


sim_dict = {'scenario' : ['hydraulics', 'hydrology', 'both'],
            "p" : [2, 4, 6, 8, 10],           
            "tr" : [10],
            "Ks_v" : [2, 5, 10, 15],
            "So" : [0.01],
            "r" : [2, 5, 10],
            "alpha_v" : [0.2, 0.5, 0.8]
            }

sim_vars = sim_dict.keys()

sim_list = [dict(list(zip(sim_vars, prod))) for prod in
                  it.product(*(sim_dict[var_name] 
                               for var_name in sim_vars))]

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
    if sim_name in os.listdir(out_dir):
        print ("skipping" , sim_name)
        return

    # shutil.rmtree(sim_dir, ignore_errors=True)  
      
    os.mkdir (sim_dir) if not os.path.isdir(sim_dir) else 0
    
    if not os.path.isdir((sim_dir+ "/Outputs")):
        os.mkdir (sim_dir + "/Outputs") 
    if not os.path.isdir((sim_dir+ "/Inputs")):
        os.mkdir (sim_dir + "/Inputs") 
    
    # update parameter dictionary for this simulation case
    case = default.copy()
    case.update(sdict) 

    D = case["D"]
    phi_veg =  case["phi_veg"]
    So =  case["So"]
    
    if case["scenario"] == 'hydraulics': 
        
        if case["Ks_v"] > 5:

            print ("skipping  " + sim_name)
            return

        case["Ks_v"] = case["Ks_b"]        

    elif case["scenario"] == 'hydrology': 
                
        case["alpha_v"] = case["alpha_b"] 
        
        # only run once ...
        if case["alpha_v"] > 0.1: 
            print ("skipping  " + sim_name)
            return
            

    case['T'] = (case['tr']+ case['t_rec'])*60  
    case['L'] = np.round(2*case['r'] + 20.4, 1)
    case['l'] = np.round(2*case['r'] + 30.4, 1)

    p0 = Params(".", case, overwrite = 1,  T = case['T']) 

    if "dt" in p0.__dict__:
        case["nbtimes"] =  p0.T/p0.dt + 1  
    else:
        case["dt"] =  p0.T/(p0.nbtimes + 1)  


    # update parameter dictionary for this simulation case
    
    x,y,z = write_planar_topo(p0, sim_dir)

    if p0.veg_type == "uniform":
        veg = write_veg_uniform(p0, sim_dir)
    
    elif p0.veg_type == "blob":
        veg = write_veg_blob(p0, sim_dir)

    elif p0.veg_type == "dot":
        veg = write_veg_dot(p0, sim_dir)        

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

    case["Nxcell"] = p0.Nxcell
    case["Nycell"] = p0.Nycell       
      
    case["imaxcoef"] = p0.Ks_v/3.6e5 # maximum infiltration rate
    write_input(case, case_dir = ".", sim_dir = sim_dir)        

    start = time.time()
    exec_name = "/Users/octaviacrompton/Dropbox/FullCSWOF/bin/FullSWOF_2D" # ../../bin/FullSWOF_2D    
    cmdout = os.popen("cd {0}; {1} ".format(sim_dir, exec_name)).read()
    # print (cmdout)
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

    #txt = input("Delete previous bank? y/n ")
    txt = "y"
    if txt != "y":
        print ("Adding to existing files")
    else:
        # remove all prev simulations
        print ("Removing existing files")
        shutil.rmtree(out_dir, ignore_errors=True) 
        os.mkdir (out_dir) if not os.path.isdir(out_dir) else 0
    
    shutil.copy(__file__, out_dir)

    os.mkdir(out_dir + "/code") if not os.path.isdir(out_dir + "/code") else 0
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