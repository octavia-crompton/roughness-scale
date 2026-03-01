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

meta = pd.read_csv(
    '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/meta_nstep20.csv')
project_dir = "/Users/octaviacrompton/Dropbox/FullCSWOF/Tests"

default = {"L" : 1, 
           "l" : 12., 
           "dx" : 1,                                
           "Ks_b" : 0.1,                       
           "p" : 5,
           "tr" : 5,           
           "dt" : 30,
           "flux" : 1,
           "veg_type" : "uniform",
           "seed" : 1,                    
           "fV" : 1,                               
           "offset" : 0,                                 
           "scheme" : 'dw',
           "fric" : 2,
           "t_rec" :  5,
           'B_bc_init' : 2,
           'imax_scale' : 2,
           'Psicoef' : 1,
            'inf' : 3
          }


def call_SWOF(args ):
    """
    Wrapper function executes fullSWOF
    """       
    sdict, key = args

    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])
        
    out_name = "calibrate_dw__7/{0}".format(key)
    out_dir = os.path.join(project_dir, out_name)        

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

    # alpha_v = case['K']# *1e-6/8
    alpha_v = case['dw']# *1e-6/8

    case['alpha_v'] = alpha_v
    case['alpha_b'] = alpha_v
    
    rain_file = '/Users/octaviacrompton/Google_Drive_quatratavia/data/polyakov2017rainfall/rain/{0}.csv'.format(key)
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


    for label in meta.label:

        row = meta.query("label == {0}".format(label)).iloc[0]
        plot = row.Plot
        ID = row.ID
        key = '{0}-{1:.0f}-{2:.0f}'.format(ID, plot, label)      
        
        out_name = key
        out_dir = os.path.join(project_dir , "calibrate_dw__7", out_name)
        folder =  os.path.join(project_dir , "calibrate_dw__7")

        if out_name in os.listdir(folder):
            print ("skipping" , out_dir)
            continue

        print(out_dir)

        os.mkdir (out_dir) if not os.path.isdir(out_dir) else 0
        
        shutil.copy(__file__, out_dir)

        os.mkdir(out_dir + "/code")
        shutil.copy("write_SWOF.py", out_dir  + "/code")
        shutil.copy("read_SWOF.py", out_dir  + "/code")

        # print("starting simulations")
        # call_SWOF(sim_list[0])    

        dw = row.f_avg
        dw = dw*0.7
        dw = np.round(np.linspace(dw*0.25**2, dw*1.5**2, 15), 2)
        dw = list(dw) 
        

        Ks_v = eval(row['Ks_S_list']).copy()
        Ks_v.append(row['Ks'])
        Ks_v = Ks_v + eval(row['Ks_R_list'])

        H = eval(row['H_S_list']) + [0]*4
        alpha_R = [10000]*4 + eval(row['alpha_R_list']).copy()

        Ks_v= [np.round(max(k, 0.1), 2).astype(float) for k in Ks_v]
        H = [np.round(h, 2).astype(float) for h in H]
        alpha_R = [np.round(k, 2).astype(float) for k in alpha_R]

        try:
           So = row['Slope']
        except:
            So = 0.1

        So = np.round(So, 3)

        sim_dict = {
                    "dw" : dw,        
                    "Ks_v" : Ks_v, 
                    "H" : H,
                    "alpha_R" : alpha_R
                    }

        sim_vars = [s for s in sim_dict.keys() if s != 'dw'] 


        sim_list0 = [(dict(list(zip(sim_vars, prod)))) for prod in
           zip(*(sim_dict[var_name] for var_name in sim_vars))]

        sim_list= []
        for i, dw in enumerate(sim_dict['dw']):  
            sim_list = sim_list + [dict(item, **{'dw':dw, 'So' : So}) for item in sim_list0]
        sim_list = [(s, key) for s in sim_list]

        print (key)

        path = os.path.join(out_dir, "sim_list.json")
    
        with open(path, 'w') as fout:
            json.dump(sim_list , fout)

        path = os.path.join(out_dir, "default.json")
        with open(path, 'w') as fout:
            json.dump(default , fout)

        start = time.time()
        from multiprocessing import Pool
        pool = Pool()
        res = pool.map(call_SWOF, sim_list)
        pool.close()
        dt = time.time() - start                  
        print (dt/60)

