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
# AUTOMATE THIS IN LOOP...
project_dir = "/Users/octaviacrompton/Dropbox/FullCSWOF/Tests"
label = 6
key = 'Ab1-{0}'.format(label)

out_name = "calibrate_optim_{0}".format(key)
out_dir = os.path.join(project_dir, out_name)


default = {"L" : 2, 
           "l" : 6., 
           "dx" : 2,                                
           "alpha_b" : 0.03,
           "alpha_v" : 0.1,
           "Ks_b" : 0.1,                       
           "p" : 5,
           "tr" : 5,           
           "dt" : 30,
           "flux" : 1,
           "veg_type" : "uniform",
           "seed" : 1,                    
           "fV" : 1,                               
           "offset" : 0,                                 
           "scheme" : 'DW',
           "fric" : 2,
           "t_rec" :  5,
           'B_bc_init' : 2,           
            'Psicoef' : 0.1,
            'imax_scale' : 2
          }

meta = pd.read_csv(
    '/Users/octaviacrompton/Dropbox/data/polyakov2017rainfall/rain/meta.csv')
meta = meta.query("label =={0}".format(label)).iloc[0]

try:
    default['So'] = meta['Slope']
except:
    default['So'] = 0.1


def call_SWOF(x0, default = default):
    """
    Wrapper function executes fullSWOF
    """       
    # path = os.path.join(out_dir, "sim_list.json")
	
    # with open(path, 'w') as fout:
    #     json.dump(sim_list , fout)

    
    sdict = {"alpha_v" : x0[0], 
            "Ks_v" : x0[1], 
            "dthetacoef" : x0[2]}

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
    
    rain_file = '/Users/octaviacrompton/Dropbox/data/polyakov2017rainfall/rain/{0}.csv'.format(key)
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

    res = read_summary(sdict, case, out_dir)
    sim = pd.Series(res.__dict__)

    from scipy import interpolate

    f = int(sim["dt"]/np.diff(sim["time"])[0])    
    diff_vol_bound = np.diff(sim.Vol_bound_tot[::f], 
        prepend = 0)/np.diff(sim['time'][::f], prepend = 0)

    qL = diff_vol_bound/sim.l/sim.L*3.6e5
    sim.hydro = qL[:len(sim.t)]

    hydro = np.insert(sim.hydro, -1, 0)
    t = np.insert(sim.t, -1, sim.t[-1] + 600)
    f = interpolate.interp1d(t/60, hydro)

    df = pd.read_csv('/Users/octaviacrompton/Dropbox/data/polyakov2017rainfall/rain/{0}.csv'.format(key))    
    runoff = f(df.Time)   # use interpolation function returned by `interp1d`

    # plt.plot(sim.t/60, sim.hydro, 'o')
    Q_rmse = np.mean((df.Runoff/10 - runoff)**2)

    print (Q_rmse)
    # import matplotlib.pylab as plt
    # plt.plot(sim.t/60, sim.hydro, '--')
    # plt.plot(df.Time, df.Runoff/10)    
    # plt.show()
    # need to process output here – read 

    print("finishing: " + sim_name)
    return  Q_rmse


def process_cmdout(cmdout):
  
    return [l for l in cmdout.split("\n") if 
            ("%] done" not in l) and (len(l)>0)]

if __name__ == "__main__":

    # txt = input("Delete previous bank? y/n ")
    txt = "y"
    if txt != "y":
        print ("Adding to existing files")
    else:
        print ("Removing existing files")
        shutil.rmtree(out_dir, ignore_errors=True) 
        os.mkdir (out_dir) if not os.path.isdir(out_dir) else 0
    
    shutil.copy(__file__, out_dir)

    os.mkdir(out_dir + "/code")
    shutil.copy("write_SWOF.py", out_dir  + "/code")
    shutil.copy("read_SWOF.py", out_dir  + "/code")



    #dw_ff = np.round(np.linspace(meta.dw_min - 20, meta.dw_min + 20, 5 ))
    #dw_ff = list(dw_ff)

    # Ks_v = meta['Ks']/10
    # Ks_v  = np.arange(Ks_v - 2, Ks_v + 0.5, 0.5)
    # Ks_v = list(np.round(Ks_v[Ks_v >0], 2))
    

    sdict = {}
    
    guess_fV = np.round(meta.dw_min*np.random.rand(), 2)
    guess_Ks_v = np.round(meta['Ks']/10*np.random.rand(), 2)
    guess_dtheta = np.round(0.2*np.random.rand(), 2)

    ## GUESS A LIST HERE....

    sdict = {"alpha_v" : guess_fV, 
            "Ks_v" : guess_Ks_v, 
            "dthetacoef" : guess_dtheta}

    # print("starting simulations")
    # call_SWOF(sdict)
    from scipy.optimize import minimize
    x0 = [sdict['alpha_v'], sdict['Ks_v'], sdict['dthetacoef']]

    # CALL SWOF...
    # GET RESULTS AND FIGURE OUT HOW TO PAIR 
    # COMPUTE OTHER EFFICIENCY METRICS FOR COMPARISON...
    # LOOK AT CRUSTED SOILS?
    # GLENN: IN DARK AGES IN TERMS OF PRESCRIBING SURFACE ROUGHNESS FOR OVERLAND FLOW...
    # A PRIORI, WITHOUT CALIBRATING.. MAKING SENSE OF THIS DATASET FROM WALNUT GULCH IMPLICATIONS BEYOND DRYLANDS, BUT 
    # SINCE SUCH A BROAD RANGE OF  SCENARIOS...

    # CROSS PARAMETERS... 

    # CUT MINIMIZE
    # RUN SAME NUMBER SIMULATIONS AGAIN...
    result = minimize(call_SWOF, x0, args = (default),
                          method='nelder-mead',
                          options= {"maxiter": 10,
                                    "return_all" : False, 
                                    # "fatol" : fatol,
                                     # "xatol" : 1, 
                                    "adaptive" : True})     

    print (result)
    # start = time.time()
    # from multiprocessing import Pool
    # pool = Pool()
    # res = pool.map(call_SWOF, sim_list)
    # pool.close()
    # dt = time.time() - start                  
    # print (dt/60)

