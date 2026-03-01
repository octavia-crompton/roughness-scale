"""
Runs fullSWOF using a template located in folder `case_name`,
and saves output to a folder `out_dir`

Domain: 2D, with gaussian-filter generated blobs 
- Vary roughness contrast

Tests skill in setting Cd - between James and James_fix

"""
import os
import sys
import shutil
import itertools as it
import contextlib
import json
from multiprocessing import Pool

my_modules = ["write_SWOF", "read_SWOF"]
for mod in my_modules:
    if mod in sys.modules: 
        del sys.modules[mod]

from write_SWOF import *
from read_SWOF import read_summary
from resistance_functions import *

# project directory
# key = 'ER2-2-108'
# key = 'TA-2-499' 
# key = 'ER4S-1-242'

# 491, 494
label = 491
meta = pd.read_csv(
    '/Users/octaviacrompton/Dropbox/data/polyakov2017rainfall/meta_missed.csv')
meta = meta.query("label =={0}".format(label)).iloc[0]
plot = meta.Plot
ID = meta.ID
key = '{0}-{1:.0f}-{2}'.format(ID, plot, label)  

project_dir = "/Users/octaviacrompton/Dropbox/FullCSWOF/Tests"

out_name = "calibrateU_{0}_L".format(key)
out_dir = os.path.join(project_dir, out_name)

default = {"L" : 0.5, 
           "l" : 6., 
           "dx" : 0.5,                                
           "Ks_b" : 0.1,                       
           "p" : 5,
           "tr" : 5,           
           "dt" : 30,
           "flux" : 1,
           "veg_type" : "uniform",
           "seed" : 1,                    
           "fV" : 1,                               
           "offset" : 0,                                 
           "scheme" : 'laminar',
           "fric" : 3,
           "t_rec" :  5,
           'B_bc_init' : 2,
           'imax_scale' : 2, 
           'Psicoef' : 0.1
          }

K = meta.K
K = list(np.round(np.linspace(K/4, K*2, 8)))

# Ks_min = min(meta['Ks_cal'],  meta['Ks'])
# Ks_max = max(meta['Ks_cal'],  meta['Ks'])
# Ks_v  = np.arange(Ks_min, Ks_max + 0.5, 0.5)
# Ks_v = list(np.round(Ks_v[Ks_v >0], 2))


Ks_v = meta['Ks']
Ks_v  = np.arange(Ks_v - 1, Ks_v + 0.5, 0.5)
Ks_v = list(np.round(Ks_v[Ks_v >0], 2))

print (len(Ks_v))

try:
    default['So'] = meta['Slope']
except:
    default['So'] = 0.1

sim_dict = {"K" : K, 
            "Ks_v" : Ks_v,
            "dthetacoef" : [0, 0.05, 0.1, 0.2]
            }

# psi = average suction across the wetting fron
# theta = moisture deficit            
sim_vars = sim_dict.keys()

sim_list = [dict(list(zip(sim_vars, prod))) for prod in
                  it.product(*(sim_dict[var_name] 
                               for var_name in sim_vars))]

def call_SWOF(sdict):
    """
    Wrapper function executes fullSWOF
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
    case = default.copy()
    case.update(sdict) 
    
    So =  case["So"]
    
    alpha_v = case['K']*1e-6/8
    case['alpha_v'] = alpha_v
    case['alpha_b'] = alpha_v
    
    rain_file = '/Users/octaviacrompton/Dropbox/data/polyakov2017rainfall/dry/{0}.csv'.format(key)
    precip = write_rain_csv(rain_file, sim_dir)
    
    case['T'] = (precip.Time.max()+ 10)*60
    
    p0 = Params(".", case, overwrite = 1, T = case['T']) 
    
    # update parameter dictionary for this simulation case
    
    x, y, z = write_planar_topo(p0, sim_dir)

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
    
    write_huv(p0, sim_dir)
    write_BC_closed(p0, sim_dir)    
    write_BC_stop(p0, sim_dir)


    if "dt" in p0.__dict__:
        case["nbtimes"] =  p0.T/p0.dt + 1  
    else:
        case["dt"] =  p0.T/(p0.nbtimes + 1)    

    case["Nxcell"] = p0.Nxcell
    case["Nycell"] = p0.Nycell       
    case["T"] = p0.T      
    case["imaxcoef"] = p0.Ks_v*p0.imax_scale/3.6e5 # maximum infiltration rate
    write_input(case, case_dir = ".", sim_dir = sim_dir)        

    start = time.time()
    exec_name = "/Users/octaviacrompton/FullCSWOF/bin/FullSWOF_2D" # ../../bin/FullSWOF_2D    
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

    # print("starting simulations")
    # call_SWOF(sim_list[0])    

    start = time.time()
    from multiprocessing import Pool
    pool = Pool()
    res = pool.map(call_SWOF, sim_list)
    pool.close()
    dt = time.time() - start                  
    print (dt/60)

