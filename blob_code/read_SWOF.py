from scipy.stats.mstats import gmean
import time
import os
import numpy as np
import re
import pandas as pd
from write_SWOF import Params
from scipy.ndimage import gaussian_filter
from shapely import geometry
import json
import matplotlib.pylab as plt
from scipy import interpolate
from topo import reshape_topo

flatten = lambda l: [item for sublist in l for item in sublist]

name_format_r = {"Infiltrated volume" : "I_vol", 
                 "Stream volume" : "h_vol",
                 'Complete volume (Inf+Stream)' : "I_h_vol",
                 'Outflow volume at the boundaries' : "outflow_vol",
                 'Time of the computation (seconds)' : "t_comp", 
                 'CPU time (clock ticks)' : "CPU_time",
                 'Number of iterations in the algorithm' : "n_iter",
                 'Froude number' : "Fr",
                 "Volume of the rain" : "rain_vol"
                }

def move_to_output(name):

    path =  '/Users/octaviacrompton/Dropbox/FullCSWOF/Tests'
    subdir= os.path.join(path, name, "model_output")
    if os.path.exists(subdir) == False:
        os.mkdir(subdir)

    files=   os.path.join(path, name, "tr-*")
    os.popen("mv {0} {1}".format(files, subdir))

    files=   os.path.join(path, name, "p-*")
    os.popen("mv {0} {1}".format(files, subdir))

    files=   os.path.join(path, name, "scenario-*")
    os.popen("mv {0} {1}".format(files, subdir))

    files=   os.path.join(path, name, "alpha_b-*")
    os.popen("mv {0} {1}".format(files, subdir))
    


def get_name(sdict):
    
    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])
    return sim_name
      

def get_summary(sim_list, default, out_dir, subfolder = '', append = '', verbose = 0 ):
    """
    read simulations
    """    
    summary = {}

    bad_list = []
    for sdict in sim_list:
        if verbose:
            print (sdict)
        try:
            res = read_summary(sdict, default, os.path.join(out_dir, subfolder))
            summary[get_name(sdict)] = res.__dict__
            print(get_name(sdict))
        except:
            bad_list.append(sdict)

    summary = pd.DataFrame(summary).T

    # fraction leaving the domain by each boundary
    # summary["R_frac"] = summary.right/summary.rain_vol
    # summary["L_frac"] = summary.left/summary.rain_vol
    # summary["T_frac"] = summary.top/summary.rain_vol
    # summary["B_frac"] = summary.bottom/summary.rain_vol

    summary.to_pickle(out_dir + "/summary{0}.pkl".format(append))

    path = os.path.join(out_dir, "bad_list{0}.json".format(append))
    with open(path, 'w') as fout:
        json.dump(bad_list , fout)
        
    return summary, bad_list


def load_summary(sim_list, default, out_dir, append = ''):
    """
    read simulations
    """    
    summary = pd.read_pickle(out_dir + "/summary{0}.pkl".format(append))

    path = os.path.join(out_dir, "bad_list{0}.json".format(append))
    try: 
        with open(path, 'r') as fout:
            bad_list = json.load(fout)
    except:
        bad_list = []
    
    return summary, bad_list

### Functions to read the inputs -----------------------------------------

def read_lines(infile):
    """
    Returns infile as a list of lines 
    """
    fin = open(infile) 
    lines = fin.readlines()
    fin.close()
    return lines

def read_rain(sim_dir, start = 1, name = 'rain.txt'):
    """
    rain = read_rain(sim_name,start = 0, name="rain.txt")
    """
    infile = os.path.join(sim_dir, "Inputs", name)

    lines = read_lines(infile)
    rain = []
    for line in lines[start:]:
        dummy = line.strip().strip("#").split(" ")
        try:
            rain.append([float(d) for d in dummy if len(d)>0])
        except ValueError:
            print( dummy)
    rain = np.array(rain)
    
    return rain

def read_array(sim_name, start = 1, name = 'topography.txt'):
    """
    """
    infile = os.path.join(input_dir(sim_name), name)

    lines = read_lines(infile)
    coords = []
    for line in lines[start:]:
        if "#" in line:
            continue
        dummy = line.strip().strip("\t").split(" ")
        try:
            coords.append([float(d) for d in dummy if len(d)>0])
        except ValueError:
            print( dummy)
    coords = np.array(coords)
    
    return coords

def read_topo_lazy(sim_dir):
    """
    clumsy code to read topography file without
     knowing file separation (\t or space)
    """
    try:
        topo = read_topo(sim_dir, sep=  "\t")
    except ValueError:
        topo = read_topo(sim_dir,  sep=  " ")
    return topo

def read_topo(sim_dir, start = 0, sep = " "):
    """
    reads topography file ending in either .txt or .dat
    """
    input_dir =  os.path.join(sim_dir, "Inputs")
    if  "topography.txt" in os.listdir(input_dir):

        infile = os.path.join(sim_dir, "Inputs", "topography.txt")
    else:
        infile = os.path.join(sim_dir, "Inputs", "topography.dat")
        

    lines = read_lines(infile)
    coords = []
    for line in lines[start:]:
        if "#" in line:
            continue

        dummy = line.strip().split(sep)

        dummy = [float(d) for d in dummy if d not in ['', '\t']]
        if len(dummy) > 0:
            coords.append(dummy)

    coords = np.array(coords)    
    return coords    


### Functions to call SWOF  ----------------------------------------------

def process_cmdout(cmdout):

    return [l for l in cmdout.split("\n") if 
            ("%] done" not in l) and (len(l)>0)]

def popen_SWOF(sim_dir):
    """
    """
    start = time.time()
    exec_name = "/Users/octaviacrompton/FullSWOF/bin/FullSWOF_2D" # ../../bin/FullSWOF_2D

    
    cmdout = os.popen("cd {0}; {1} ".format(sim_dir, exec_name)).read()

    end = time.time()
    return process_cmdout(cmdout),  (end-start)


### Functions to read the output ------------------------------------------


def read_summary(sdict, default, out_dir):
    """
    Added 2-d fields:
        norm_U : overland flow velocity (u**2 + v**2)**0.5 at t=tr
        infl_2d : cumulative infiltration after storm, units = m 
            (multiply by dx**2 to get volume)  
    """ 
    sim_name = ','.join(['-'.join([key, str(sdict[key])])
                               for key in list(sdict.keys())])

    sim_dir = os.path.join(out_dir, sim_name)
    
    results = read_results(sim_dir)

    default.update(sdict)

    p0 = Params(sim_dir + "/Inputs", default, mode = '')
    
    p0.results = results

    for key in name_format_r:
        setattr(p0, name_format_r[key], results[key])    

    # alpha = read_alpha(get_name(sdict), p0, out_dir)
    # p0.alpha = alpha
    # p0.alpha_avg = alpha.mean()

    veg = read_veg(get_name(sdict), out_dir, p0.Nxcell, p0.Nycell)

    p0.veg = veg
    p0.veg_avg = veg.mean()    

    Ks = read_Ks(get_name(sdict), p0, out_dir)
    p0.Ks = Ks

    if p0.veg_type != 'uniform' and veg.std() ==0:
        # This case reads vegetation from the Ks file
        veg = (Ks > Ks.min()).astype(float)
        # print (veg.mean())
        p0.veg = veg
        p0.veg_avg = veg.mean()    

    check_vol = read_check_vol(sim_name, out_dir)
    
    evol = process_evolution(sim_name, p0, out_dir)
    
    for key in evol.keys():
        setattr(p0, key, evol[key])     

    for key in check_vol.keys():
        setattr(p0, key, np.array(check_vol[key]))

    # try:
    #     BT = read_BT(sim_dir, p0)
    #     setattr(p0, "BT", BT)
    # except:
    #     pass

    if p0.rain_vol == 0:
        setattr(p0, "rain_vol", evol["hc"][0].sum()/p0.dx**2) 

    p0.IF = p0.I_vol/p0.rain_vol
    p0.infl_frac = p0.I_h_vol/p0.rain_vol
    p0.t_rain = p0.tr*60
    p0.i_tr = int(p0.t_rain/p0.dt)

    p0.Vol_bound_tot = np.array(check_vol["Vol_bound_tot"])
    p0.time = np.array(check_vol["time"])

    topo = read_topo_lazy(sim_dir)
    xc, yc, zc = reshape_topo(topo, p0)
    p0.xc = xc
    p0.yc = yc
    p0.zc = zc
    boundaries = read_boundary(p0, sim_name, out_dir)

    # cumulative 
    for key in [ 'left', 'right', 'bottom', 'top', 'outflow_m3']: 
        setattr(p0, key + "_c", boundaries[key].iloc[-1])     

    return  p0


def read_results(sim_dir):
    """
    reads "results.dat"
    """
    infile = os.path.join(sim_dir, 'Outputs', 'results.dat')

    lines = read_lines(infile)
    lines = [l.strip() for l in lines if (("#" not in l) & (len(l)>1))]
    results = []
    lines = [t.split(":") for t in lines]
    for i, l in enumerate(lines):
        lines[i][1] = float(lines[i][1])
        
    # Rename : mydict[new_key] = mydict.pop(old_key)
    return dict(lines)

def percent_diff(s, match, p0):    
    """
    percent difference in hydrographs between two simulations
    """
    tmin = max(s.time[0], match.time[0])
    tmax = min(s.time[-1], match.time[-1] )
    t = np.arange(tmin, tmax)
    f = interpolate.interp1d(s.time, s.Vol_bound_tot)
    q_swof = f(t)       
    f = interpolate.interp1d(match.time, match.Vol_bound_tot)    
    q_match = f(t)      
    hydro_percent_diff = np.zeros_like(q_swof)
    
    inds =  np.where(q_swof > p0.p_mps*0.05)[0]
    hydro_percent_diff[inds] = (q_swof[inds] - q_match[inds])/q_match[inds]*100.
    hydro_percent_diff = hydro_percent_diff.mean()
    
    return (hydro_percent_diff)

def process_evolution(key, p0, out_dir):
    """

    """
    df =  read_evolution(key, out_dir)    

    t =  np.array(df.t.unique())

    evol = {}

    h_1d = np.array(df['h']).reshape(
            len(t), p0.Nxcell, p0.Nycell).mean(1).mean(1)
    
    try:
        end = np.where(h_1d > 1e-7)[0][-1]+5
    except IndexError:
        end = 0

    for key in ["h", "u", "v"]:
        evol[key+ "c"] = np.array(df[key]).reshape(
            len(t), p0.Nxcell, p0.Nycell)#[:end]

    for key in ["I"]:
        evol[key] = np.array(df[key]).reshape(
            len(t), p0.Nxcell, p0.Nycell)#[:end]

    evol["t"] = t#[:end]
    evol["infl_2d"] = evol["I"][-1]

    try:
        evol["h_max"] = evol["hc"][p0.i_tr]
        evol["q_max"] = evol["vc"][p0.i_tr]*evol["hc"][p0.i_tr]
    
    except:        
        evol["h_max"] = evol["hc"][1]
        evol["q_max"] = evol["vc"][1]*evol["hc"][1]
        evol["i_tr"] = 1

    
    return evol


def read_veg(sim_name, out_dir, Nxcell, Nycell, start = 1, name = 'veg.txt'):
        """
        Parameters:
        ----------
        sim_name : unique simulation key

        out_dir : simulation directory
        """
        infile = os.path.join(out_dir, sim_name, "Inputs", name)

        lines = read_lines(infile)
        coords = []
        for line in lines[start:]:
            if "#" in line:
                continue

            dummy = line.strip().strip("\t").split(" ")
            try:
                coords.append([float(d) for d in dummy if len(d)>0])
            except ValueError:
                print( dummy)
        coords = np.array(coords)
        
        return coords[:, 2].reshape((Nxcell, Nycell))


def read_Ks(sim_name, p0, out_dir, start = 1, name = 'Ks.txt'):

    infile = os.path.join(out_dir, sim_name, "Inputs", name)

    lines = read_lines(infile)
    coords = []
    for line in lines[start:]:
        if "#" in line:
            continue

        dummy = line.strip().strip("\t").split(" ")
        try:
            coords.append([float(d) for d in dummy if len(d)>0])
        except ValueError:
            print( dummy)
    coords = np.array(coords)
    
    return coords[:, 2].reshape((p0.Nxcell, p0.Nycell))


def read_evolution(sim_name, out_dir):
    """
    Read huz_evolution.dat from outputs folder
    """
    sim_dir = os.path.join(out_dir, sim_name)
    
    infile = os.path.join(sim_dir, "Outputs", 'huz_evolution.dat')

    lines = read_lines(infile)

    lines = [l.strip() for l in lines if  (len(l)>1)]

    dummy = lines[5:]
    huv = []
    for i, l in enumerate(dummy):
        if "time" in l:    
            t = (re.findall(r"[-+]?\d*\.\d+|\d+", l))[0]
        else:
            l = [col.strip() for col in l.split("\t")]
            if len(l) > 1:
                l.append(t)
            huv.append(l)
    huv = np.array(huv, dtype = float)
    
    try:
        colnames= ['x_cc', 'y_cc', 'h', 'u', 'v', 'z', 'I',  't']   
        df = pd.DataFrame(huv, columns=colnames)
    except: 
        colnames = lines[3].strip("#").strip()
        colnames = [col.strip() for col in colnames.split("\t")]
        print (huv)
        df = pd.DataFrame(huv, columns=colnames)

    return df.astype(float)


def read_huz(case_dir, name = 'huz_final.dat'):
    """
    Reads huz files

    Parameters:
    ----------
    case_dir
    
    name: specify output file
        huz_final.dat or huz_initial.dat)
    """
    infile = os.path.join(case_dir, "Outputs", name)

    lines = read_lines(infile)

    if name == 'huz_final.dat':
        colnames = lines[4].strip("#").strip()
        colnames = [col.strip() for col in colnames.split("\t")]  

        colnames = ['x_cc', 'y_cc', 'h', 'u', 'v', 'z']
    elif name == 'huz_initial.dat':   
        colnames = lines[3].strip("#").strip()
        colnames = [col.strip() for col in colnames.split("\t")]

        colnames = ['x_cc', 'y_cc', 'h', 'u', 'v', 'h+z', 'z'] 

    huz = []
    for i, l in enumerate(lines):    
        if (len(l) > 1) and ("#" not in l):
            l = [col.strip() for col in l.split("\t")]
            huz.append(l)

    return pd.DataFrame(huz, columns=colnames).astype(float)

############### Debug outputs #########################

def read_boundary(p0, key, out_dir, name = 'boundaries_flux.dat'):
    """
    NOTE: fluxes are in m^3
    """   

    infile = os.path.join(out_dir, key, "Outputs", "boundaries_flux.dat")

    lines = read_lines(infile)

    colnames  = lines[4].strip("# ").strip("\n").split("\t")    

    fluxes = []
    for line in lines:
        if "#" in line:
            continue

        dum = line.strip("\n")
        try:
            dum = [float(d) for d in dum.split("\t")]
            fluxes.append(dum)
        except ValueError:
            print( dum)
            
    fluxes = np.array(fluxes)
    
    boundary_names = {'time': 't_boundary',
        'left cumulative flux (m^2)': 'left',
        'right cumulative flux (m^2)': 'right',
        'bottom cumulative flux (m^2)': 'bottom',
        'top cumulative flux (m^2)' : 'top'
        }

    boundary = pd.DataFrame(fluxes, columns=colnames).astype(float)
    
    b = boundary.rename(columns = boundary_names)

    b["left"] = b["left"]*p0.dx
    b["bottom"] = b["bottom"]*p0.dx
    b["right"] = b["right"]*p0.dx
    b["top"] = b["top"]*p0.dx

    b["outflow_m3"] = - b.left - b.bottom + b.right + b.top
    
    return b


def read_side_boundary(sim_dir, name = "flux_boundaries_LR.dat" ):
    """
    Units are m^2/s
    """
    infile = os.path.join(sim_dir, "Outputs", name)

    lines = read_lines(infile)

    colnames  = lines[4].strip("# ").strip("\n").split("\t")    
    colnames = [c.strip() for c in colnames]
    fluxes = []
    for l in lines:
        if (len(l.strip("\n")) > 1) and ("#" not in l):
            dum = [float(d) for d in l.strip("\n").split("\t")]

            fluxes.append(dum)

    fluxes = np.array(fluxes)

    return pd.DataFrame(fluxes, columns=colnames).astype(float)


def read_check_vol(key, out_dir):
    """
    Read check_vol.dat in Outputs folder

    check_vol.dat contains the cumulated volumes [m3] at each time step. 
    
    Columns: 
         overland flow volume (Vol_of_tot) 
         the infiltrated volume (Vol_inf_tot), 
         the rain volume (Vol_rain_tot) and 
         the balance of the input and output volumes at the boundaries (Vol_bound_tot).

    Usage:  
        read_check_vol(sim.name, out_dir)
    
    """
    infile = os.path.join(out_dir, key, "Outputs", "check_vol.dat")
    lines = read_lines(infile)

    colnames = lines[3].strip("#").strip()
    colnames = [col.strip() for col in colnames.split("\t")]  

    check_vol = []
    for i, l in enumerate(lines):    
        if (len(l) > 1) and ("#" not in l):
            l = [col.strip() for col in l.split("\t")]
            check_vol.append(l)

    return pd.DataFrame(check_vol, columns=colnames).astype(float)


#### Landscape descriptors

from scipy import ndimage 

def get_edges_b(image):
    """
    edges that are bare soil areas
    """
    padded = np.pad(image, [(1, 1), (1, 1)], mode='edge')
    edges_b = (padded - ndimage.morphology.binary_dilation(padded))[1:-1, 1:-1]
    return np.abs(edges_b)

def get_edges_v(image):
    """
    edges that are vegetated areas
    """
    padded = np.pad(image, [(1, 1), (1, 1)], mode='edge')
    edges_v = (padded - ndimage.morphology.binary_erosion(padded))[1:-1, 1:-1]
    return edges_v



### More processing of the SVE output
def get_alpha(sim):
    """
    get the simulation alpha field
    """
    alpha = sim.veg.copy().astype(float)
    alpha[sim.veg < 1] = sim.alpha_b
    alpha[sim.veg == 1] = sim.alpha_v
    
    return alpha


def forward_append(temp, sim):
  
    while len(temp) < len(sim.t):
        temp = np.append( [0], temp)  
        
    while len(temp) > len(sim.t):
        temp = temp[1:]        
    
    return temp

def backward_append(temp, sim):
    
    while len(temp) < len(sim.t):
        temp = np.append(  temp, [0])  

    while len(temp) > len(sim.t):
        temp = temp[:-1]
        
    return temp


def add_cmo(summary):
    cmos = []
    for key in summary.index:
        sim = summary.loc[key]
        
        cmo = ndimage.measurements.center_of_mass(sim.veg)[1]
        cmos.append(cmo)

    summary["cmo"] = cmos
    return summary


def read_BT(sim_dir, p0):
    """
    To compare output to the hydrograph:
        
        BT_sum = BT.groupby("time").sum().reset_index()
        plt.plot(BT_sum["time"], BT_sum["ctop"])
        plt.plot(sim["time"], sim.Vol_bound_tot)

    To convert between units:

        BT2 = BT.query("x_cc == 1.25")
    
        plt.plot(BT2["time"][1:], np.diff(BT2.ctop)/sim["dt"]/sim.dx)
        plt.plot(BT2["time"], BT2.top_m2s)


    """
    BT = read_side_boundary(sim_dir, name = "flux_boundaries_BT.dat" )
    dt = np.diff(BT.time.unique())[0]

    BT.columns = ['x_cc', 'bottom_m2s', 'top_m2s', 'time']

    BT = BT.groupby(["time", "x_cc"]).max().reset_index()
    try:
        f = int(p0.dt/np.diff(p0.time)[0])
    except:
        f = int(p0["dt"]/np.diff(p0.time)[0])

    BT[['cbottom', 'ctop']] = BT.groupby(["x_cc"]).cumsum()[['bottom_m2s', 'top_m2s']]*dt*p0.dx
    
    
    keep = list(set([t for t in BT["time"] if t in p0.t]))

    keep.sort()
    inds = [i for i in range(len(BT)) if BT.time[i] in keep] 
    if len(BT.query("x_cc == {0}".format(np.min(BT.x_cc)))) != len(p0.hc):
        inds +=  list(BT.query("time == {0}".format(np.max(BT["time"]))).index)

    BT = BT.loc[inds]

    return BT


def insert_fld(summary, fld):
    """
    Add object field to pandas dataframe
    """
    summary[fld] = 0
    summary[fld] = summary[fld].astype(object)
    
    return summary

def instabilities(summary, threshold = 1):


    summary = insert_fld(summary,"uc_raw")
    summary = insert_fld(summary,"vc_raw")    

    for key in summary.index:
        
        sim = summary.loc[key]
        
        if not threshold:
            threshold = np.percentile(sim.vc[sim.vc > 0], 99.99)
        
        u_inst = (np.abs(sim.uc) > threshold).sum()

        if u_inst > 0:

            summary.at[key, 'u_inst'] = u_inst

            U = sim.uc.copy()
            # print (sim.name, u_inst)

            t, x, y = np.where(np.abs(sim.uc) > threshold)

            yy = y[y < sim.Nycell - 1]
            xx = x[y < sim.Nycell - 1]
            tt = t[y < sim.Nycell - 1]    
            
            if len(yy) > 0:
                U[tt, xx , yy] =   (U[tt, xx , yy - 1]  + U[tt, xx , yy + 1] )/2.
                
            yy = y[y == sim.Nycell - 1]
            xx = x[y == sim.Nycell - 1]
            tt = t[y == sim.Nycell - 1]            

            if len(yy) > 0:
                U[tt, xx , yy] =   U[tt, xx, yy - 1]

            #summary.at[key, 'uc_raw'] = sim.uc
            summary.at[key, 'uc'] = U


        v_inst = (sim.vc > threshold).sum()            

        if  v_inst> 0:

            summary.at[key, 'v_inst'] = v_inst            

            U = sim.vc.copy()

            t, x, y = np.where(sim.vc > threshold)
            xx = x[x < sim.Nxcell - 1]
            yy = y[x < sim.Nxcell - 1]
            tt = t[x < sim.Nxcell - 1]    

            if len(xx) > 0:
                U[tt, xx , yy] =   (U[tt, xx - 1, yy]  + U[tt, xx + 1, yy] )/2.

            xx = x[x >=  sim.Nxcell - 1]
            yy = y[x >= sim.Nxcell - 1]
            tt = t[x >= sim.Nxcell - 1]     
            if len(xx) > 0:
                U[tt, xx , yy] =  U[tt, xx - 1, yy] 

            # summary.at[key, 'vc_raw'] = sim.vc
            summary.at[key, 'vc'] = U


    return summary


def compute_mass_bal(summary):
    
    for key in summary.index:
        sim = summary.loc[key]

        mb = (sim.Vol_rain_tot[-1] - sim.Vol_bound_tot[-1] - sim.Vol_of_tot[-1] - sim.Vol_inf_tot[-1])
        summary.at[key, 'mass_bal'] = mb
        summary.at[key, 'mass_bal_percent'] = mb/sim.Vol_rain_tot[-1]*100    

    return summary
    
def add_Fr(summary):
    """
    add froude number
    """
    for key in summary.index:

        sim = summary.loc[key]
        Fr = sim.vc[sim.i_tr]/np.sqrt(sim.hc[sim.i_tr]*9.8)
        Fr[sim.hc[sim.i_tr] < 1e-4] = np.nan
        Fr = Fr[~np.isnan(Fr)]
        if len(Fr) == 0:
    
            summary.at[key, 'Fr'] = np.nan
            summary.at[key, 'Fr_max'] = np.nan
    
        else:
    
            if (np.nanmax(Fr)) > 1:
                print (sim.name,  np.nanmax(Fr))                     
    
            summary.at[key, 'Fr'] = np.nanmean(Fr)
            summary.at[key, 'Fr_max'] = np.percentile(Fr, 99)
    
    return summary  

def add_Re_all(summary):
    """
    over the storm duration
    """
    for key in summary.index:
        sim = summary.loc[key]
        Re = sim.vc*sim.hc/1e-6

        summary.at[key, 'Re_all'] = np.mean(Re)
        
    return summary

def add_terms(summary):
    """
    Add hydrograph, SVE terms (with RMSE measures) to summary dataframe
    
    Assumes Mannings roughness

    Derived fields  have the same shape as t,hc,uc,vc,I
    
    hydro: hydrograph, cm/hr   [nbtimes]
    d_hydro: rate of change of hydrograph, cm/hr^2  [nbtimes]

    final_qL : hydro[-1], cm/hr; tests whether simulations were cut short [float]

    I_v : infiltration in vegetation, normalized by rainfall amount
        no run-on : 0
        accounts for vegetation fraction
    
    I_v_tr  : infiltration in vegetation prior to the end of the rain, normalized by rainfall amount


    IF_v  : fraction of rainfall that infiltrates within vegetated areas 
        does not account for vegetation fraction

            
    """
    for fld in ['hydro', 'd_hydro', 'final_qL']:
        
        summary = insert_fld(summary,fld)

    for key in summary.index:

        sim = summary.loc[key]

        # f = int(sim["dt"]/np.diff(sim["time"])[0])    
        # f = max(f, 1)
        x = sim['time']
        y = sim['Vol_bound_tot']
        xnew = np.arange(0, sim['time'][-1], sim['dt'])
        ynew = np.interp(xnew, x, y)

        diff_vol_bound = np.diff(ynew, prepend = 0)/sim['dt']
        
        qL = diff_vol_bound/sim.l/sim.L*3.6e5
        hydro = qL[:len(sim.t)]
        
        diff_len =  len(sim.t) - len(hydro) 


        if diff_len > 0:
            hydro = np.append(hydro, [0] * diff_len)


        summary.at[key, "hydro"] = hydro

        d_hydro = np.diff(qL, prepend = 0)/sim['dt']*3600
        summary.at[key, "d_hydro"] = d_hydro[:len(sim.t)]
        summary.at[key, "final_qL"] = qL[-1]
        inds = sim.hc < 5e-4
        try:
            summary.at[key, 'flashy'] = sim.t[np.where(hydro > 0.99*hydro.max())[0][0]]/60 # min  
        except: 
            summary.at[key, 'flashy'] =  np.nan

        # steady =  (sim.p - sim.fV*sim.Ks_v  - qL.max())/(sim.p)*100
        # summary.at[key, 'steady'] = steady.round(3)
        
        if sim.Bbound == 5:
            
            summary.at[key, 'I_v'] = sim.infl_2d[sim.veg == 1].sum()*sim.dx**2/(
                        sim.bottom_c +  sim.rain_vol)  
    
            summary.at[key,'IF_v'] = sim.infl_2d[sim.veg == 1].sum()*sim.dx**2/(
                sim.bottom_c + sim.rain_vol) 

            summary.at[key,'I_v_tr'] =  sim.I[sim.i_tr, sim.veg == 1].mean()*100/(
                sim.bottom_c +  sim.rain_vol)  

        else:
            summary.at[key,'I_v'] = (sim.infl_2d[sim.veg == 1].mean()*100)/(sim.p*sim.tr/60) 
            
            summary.at[key,'IF_v'] = sim.infl_2d[sim.veg == 1].sum()*sim.dx**2/sim.rain_vol
        
            summary.at[key,'I_v_tr'] = sim.I[sim.i_tr,sim.veg == 1].mean()*100/(sim.p*sim.tr/60) 
      

        summary.at[key, 'C'] = 1 - sim.infl_frac

    try:
        summary['excess'] = (summary.fV*(summary.p - summary.Ks_v) +
                     (1 - summary.fV)*( summary.p -  summary.Ks_b))
        summary['excess_tr'] = summary['excess']*summary['tr']/60

    except:
        summary['excess'] = 0
        summary['excess_tr'] = 0
    
    summary['p-Ks_v'] = summary['p'] - summary['Ks_v']
    summary['I_v_rec'] =   summary['I_v'] - summary['I_v_tr']
    
    return summary 

def get_flowlines(sim):
    """
    """

    flow_lines = []
    for i in range(1, sim.Nxcell):
        starting_point = np.array([[0.05, i*0.05]])

        strm = plt.streamplot(sim.yc[0], sim.xc[:, 0], 
                              sim.vc[sim.i_tr], sim.uc[sim.i_tr], 
                              density=1, start_points=starting_point, 
                              color = 'green');


        num_pts = len(strm.lines.get_segments())
        flow_line = np.full((num_pts, 2), np.nan)
        for i in range(num_pts):
            flow_line[i,:] = strm.lines.get_segments()[i][0,:]

        flow_lines.append(flow_line)
    return flow_lines

def add_flowlines(summary):
    """
    """
    for key in summary.index:
        sim = summary.loc[key]
        flow_lines = get_flowlines(sim)

        y_pos = flatten([list((fl[:, 1]  - fl[:, 1].mean(0))) for fl in flow_lines])

        summary.at[key, 'std_flowlines'] = np.std(y_pos)

    return summary

def add_U_max(summary):
    """ 
    v (code) <-> U (notation)
    u (code) <-> V (notation)    
    """ 
    for key in summary.index:
        
        sim = summary.loc[key]
        
        summary.at[key, 'U_tr_mean'] = sim['vc'][sim.i_tr].mean()
        summary.at[key, 'V_tr_mean'] = sim['uc'][sim.i_tr].mean()
        
        summary.at[key, 'U_max'] = sim['vc'][sim.i_tr].mean()
        summary.at[key, 'U_std'] = sim['vc'][sim.i_tr].std()
        #summary.at[key, "U_std"] = sim.vc.std()    
        
        summary.at[key, "V_std"]  = sim.uc[sim.i_tr].std()
        
        summary.at[key, "std(V/U)"]  = (sim.uc/(sim.vc+ 1e-4)).std()
        #summary.at[key, "V_std"] = (sim.uc).std()
        
    return summary


def add_kinematic_terms(summary):
    """
    Add hydrograph, SVE terms (with RMSE measures) to summary dataframe
    
    *** Update and do not assume Manning's roughness ***
    U_KWE : U_KWE(h) in m/s; U predicted using KWE  [nbtimes]
    dU_KWE : U -U_KWE(h) in m/s; U  difference from kinematic [nbtimes]

    dSf_KWE : Sf - So; Sf difference from kinematic [nbtimes]

    local : dU/dt; m/s, local acceleration term
     
    ADD COMMENTS
    
    """
    for fld in [ 'U_KWE', 'dU_KWE',
                 'dSf_KWE',  'local', 'convective', 'pgrad', 'vertical']:
        
        summary = insert_fld(summary,fld)

    for key in summary.index:

        sim = summary.loc[key]

        f = int(sim["dt"]/np.diff(sim["time"])[0])    
        inds = sim.hc < 5e-4
        
        # 1/g <dU/dt>  [-]
        dt = np.diff(sim["t"])
        dt = dt.reshape([len(dt), 1,1])
        
        t1 = np.diff(sim.vc, axis = 0)/dt/9.8
        t1 = np.concatenate((t1, t1[-1:]), axis = 0)
        t1[inds] = 0
        local = t1.mean(1).mean(1)
        summary.at[key, "local"]  = local

        # < U/g dU/dx>     [-]   
        t2 = np.diff(sim.vc, axis = 2, prepend=sim.vc[:, :, :1]*0)
        t2 = (sim.vc*t2)/sim.dx/9.8
        t2[inds] = 0
        convective = t2.mean(1).mean(1)        
        summary.at[key, "convective"]  = convective

        #  <dh/dx>    [-]   
        t3  =  np.diff(sim.hc, axis = 2, prepend=sim.hc[:, :, :1]*0)/sim.dx
        t3[inds] = 0
        pgrad = t3.mean(1).mean(1)        
        summary.at[key, "pgrad"] = pgrad
        
        #  1/g <(U/h) (p - I)>  [-]          
        diff_rain_bound = np.diff(sim.Vol_rain_tot[::f])/np.diff(sim.time[::f])  # m3/s
        p = forward_append( diff_rain_bound/sim.L/sim.l, sim) # m/s
        p = p.reshape([len(p), 1, 1])

        I = np.diff(sim.I, axis = 0)/dt/sim.dx**2
        I = np.concatenate((I, I[-1:]), axis = 0)
        t4 = (sim.vc/sim.hc*(p-I))/9.8
        t4[inds] = 0
        vertical = t4.mean(1).mean(1)
        summary.at[key, "vertical"] = vertical
  
        # summary.at[key, "dU_KWE_rmse"] =  np.mean(np.sqrt(dU_KWE**2))*100
        # summary.at[key, "dSf_KWE_rmse"] =  np.nanmean(np.sqrt(dSf_KWE**2))/sim.So  
        summary.at[key, "local_rmse"] =  np.nanmean(np.sqrt(local**2))/sim.So 
        summary.at[key, "convective_rmse"] =  np.nanmean(np.sqrt(convective**2))/sim.So       
        summary.at[key, "pgrad_rmse"] =  np.nanmean(np.sqrt(pgrad**2))/sim.So           
        summary.at[key, "vertical_rmse"] =  np.nanmean(np.sqrt(vertical**2))/sim.So               

    return summary 

###################################### Code to describe hysteresis ###################################

def get_hydro_hyst_orient(sim):
    """
    Q(L) vs dQ(L)/dt hysteresis 

    hydro_hyst_orient < 0 means q_L(<h>) is larger 
        during the rising limb than the falling limb (clockwise)
    
    """
    from scipy import interpolate

    f = int(10/np.unique(np.diff(sim.time))[1].round(2))
    x = sim.Vol_of_tot[::f]

    y = np.diff(sim.Vol_bound_tot[::f], prepend = [0])

    t = sim.time[::f]
    x = x/x.max()
    y = y/y.max()

    r_inds = t <= sim.tr*60
    f_rise  = interpolate.interp1d(x[r_inds],y[r_inds])

    f_inds = t >= sim.tr*60
    f_fall = interpolate.interp1d(x[f_inds],y[f_inds])

    x_min = np.maximum(x[f_inds].min(), x[r_inds].min())
    x_max = np.minimum(x[f_inds].max(), x[r_inds].max())
    x_interp = np.arange(x_min, x_max, 0.01)
    q_rise = f_rise(x_interp)      
    q_fall = f_fall(x_interp)      


    return np.trapz( q_fall - q_rise , x_interp)



def describe_hystersis(summary):
    for key in summary.index:
        sim = summary.loc[key]

        h_1d = sim.hc.mean(1).mean(1)
        u_1d = sim.uc.mean(1).mean(1)
        q_1d = sim.qc.mean(1).mean(1)    

        # q_hyst : <h> versus <q>
        poly = geometry.Polygon(zip(h_1d/h_1d.max(), q_1d/q_1d.max()))        
        summary.at[key, "q_hyst"] = poly.area


        # q_hyst : <h> versus Q(L)
        f = int(10/np.diff(sim.time)[0])    
        diff_vol_bound = np.diff(sim.Vol_bound_tot[::f], prepend = [0])/np.diff(sim["time"][::f], prepend = [0])
        diff2_vol_bound = np.diff(diff_vol_bound, prepend = [0])
        vol_of = sim.Vol_of_tot[::f]  

        poly = geometry.Polygon(zip(vol_of/vol_of.max(), diff_vol_bound/diff_vol_bound.max()))
        summary.at[key, "hydro_hyst"] = poly.area    

        hydro_hyst_orient = get_hydro_hyst_orient(sim)
        summary.at[key, "hydro_hyst_orient"] = hydro_hyst_orient

        poly = geometry.Polygon(zip( diff_vol_bound/diff_vol_bound.max(),  diff2_vol_bound/diff2_vol_bound.max()))
        summary.at[key, "obs_hyst"] = poly.area    

        summary.at[key, "h_std"] = sim.hc[sim.i_tr].std()/sim.hc[sim.i_tr].mean()
        summary.at[key, "q_std"] = sim.qc[sim.i_tr].std()/sim.qc[sim.i_tr].mean()

        summary.at[key, "fV_avg"] = sim.veg.mean()
        summary.at[key, "cmo_tr"] = ndimage.measurements.center_of_mass(sim.hc[sim.i_tr])[1]

        summary["clockwise"] = summary["hydro_hyst_orient"] < 0      

    return summary

###################################### Code to describe vegetation patterns ###################################

def get_patchL(veg, saturate):
    """
    compute patch lengths
    """
    nrow = veg.shape[0]
    ncol = veg.shape[1]
    patch_LV = np.zeros(veg.shape, dtype = float)  # veg patch length
    upslope_B = np.zeros(veg.shape, dtype = float)  # upslope interspace patch length (paired to veg patch)
    
    for i in range(nrow):  # loop over across-slope direction first
        count = 0           
        for j in range(ncol):  

            if veg[i, j] == 1:    #  if veg patch, add 1                
                count += 1  
                if j >= (ncol -1):  # if we're at the top of the hill   
                    patch_LV[i, j-count+1:] = count  # record veg patch length                  
                                                           
            # if [i,j] is bare and the slope cell is vegetated, record.
            # each patch starts at [i,j-count] and ends at [i,j-1]
            elif veg[i,j] == 0 and veg[i, j-1] == 1:   
                if j > 0:
                    # veg patch starts at j-count and ends at j
                    patch_LV[i, j-count:j] = count
                    try:
                        # find the nearest upslope veg cell
                        Lb = np.where(veg[i,j:] == 1)[0][0]                               
                        upslope_B[i,j-count:j] = Lb
                    except IndexError:  # bare patch extends to top of hill
                        upslope_B[i,j-count:j] = ncol - j
                    count = 0 
    patch_LV[patch_LV > saturate] = saturate
    upslope_B[upslope_B > saturate] = saturate
        
    return  patch_LV, upslope_B

def get_flowlength(veg):
    """
    compute patch lengths
    """
    nrow = veg.shape[0]
    ncol = veg.shape[1]

    flowlength = np.zeros(veg.shape, dtype = float)  
    
    for i in range(nrow):  # loop over across-slope direction first
        count = 0           
        for j in range(ncol-1, 0, -1):  
            
            if veg[i, j] == 0:    #  if not a veg patch, add 1                
                flowlength[i, j] = count  # record veg patch length                  
                count += 1  

            if veg[i,j] == 0 and veg[i, j-1] == 1:               
                flowlength[i, j] = count  # record veg patch length                  
                count = 0
        
    return  flowlength

def add_lengthscales(summary):
    """
    """
    import cv2
    for fld in ['patch_LV', 'LV_dist', 'upslope_B', 'upslope_V', 'uB_dist',
                 'patch_WV', 'WV_dist', 'patch_LB', 'LB_dist', 'L_dc_dist']:
        
        summary = insert_fld(summary,fld)

    for key in summary.index:
        
        sim = summary.loc[key]
        
        summary.at[key,"fV"] = sim.veg.mean()

        patch_LV, upslope_B = get_patchL(sim.veg, 1000) 
        
        summary.at[key,"patch_LV"] = patch_LV # [patch_LV> 0].mean()
        
        summary.at[key,"LV_dist"] =  unweight_array(patch_LV)
        summary.at[key,"LV"] =  np.mean(summary.at[key,"LV_dist"])    
                
        summary.at[key,"upslope_B"] = upslope_B
        uB = unweight_array(upslope_B)
        summary.at[key,"uB_dist"] = uB
        summary.at[key,"uB"] = np.mean(uB)
 
        patch_WV, patch_right_B = get_patchL(sim.veg.T, 1000) 
        summary.at[key,"patch_WV"] = patch_WV
        WV = unweight_array(patch_WV) 
        summary.at[key,"WV_dist"] = WV
        summary.at[key,"WV_avg"] = np.mean(WV)       
        
        patch_LB, upslope_V = get_patchL(1-sim.veg, 1000)   
        L_dc = patch_LB[:, -1]
        summary.at[key,"L_dc_dist"] =  L_dc
        summary.at[key,"L_dc"] = L_dc.mean()     
        summary.at[key,"L_dc>0"] = L_dc[L_dc>0].mean()     
        summary.at[key,"patch_LB"] = patch_LB
        summary.at[key,"upslope_V"] = upslope_V

        # bare soil length, including 0s
        summary.at[key,"LB_w"] = patch_LB[patch_LB>0].mean() # weighted 
        summary.at[key,"LB_wv"] = patch_LB.mean() # weighted with vegetation
        
        # bare soil length, unweighted
        LB = unweight_array(patch_LB) 
        summary.at[key,"LB_dist"] = LB
        summary.at[key,"LB"] = np.mean(LB) # unweighted
        
        # along slope fraction that is transition from bare soil to vegetated
        bare_to_veg = (np.diff(sim.veg, 1) == -1)
        summary.at[key, "bare_to_veg"] = bare_to_veg.mean()
        
        # # along slope fraction that is transition from vegetated to bare soil
        # (np.diff(sim.veg, 1) == 1).mean()
        veg_to_bare = (np.diff(sim.veg, 1) == 1)  
        summary.at[key,"veg_to_bare"] = veg_to_bare.mean()

        # # fraction of hillslope available for lateral/transverse flow
        # (np.diff(sim.veg, 1) != 1).mean()


        edge = (sim.veg - cv2.erode(sim.veg, np.ones([1,3]), 1)) 
        summary.at[key,'edge'] =  edge.sum()/sim.veg.sum()
    

    return summary



def add_Reynolds_numbers(summary):
    for key in summary.index:
        sim = summary.loc[key]
        Re = sim.vc[sim.i_tr]*sim.hc[sim.i_tr]/1e-6

        summary.at[key, 'Re_out'] = np.mean(Re[:, -1])
        
        summary.at[key, 'Re_90'] = np.percentile(Re.ravel(), 90)

        summary.at[key, 'Re'] = np.mean(Re)    
        summary.at[key, 'Re_v'] = np.mean(Re[sim.veg == 1])            
        summary.at[key, 'Re_b'] = np.mean(Re[sim.veg == 0])                    

        Fr2 =  (sim.vc[sim.i_tr]**2/(sim.hc[sim.i_tr]*9.8))
        dw = 8/Fr2
        dw = 8*(sim.hc[sim.i_tr]*9.8)/sim.vc[sim.i_tr]**2
        dw[sim.vc[sim.i_tr] < 1e-2] = np.nan
        summary.at[key, 'dw'] = np.nanmean(dw)
        summary.at[key, 'dw_v'] = np.nanmean(dw[sim.veg == 1])
        summary.at[key, 'dw_b'] = np.nanmean(dw[sim.veg == 0])

    return summary


def unweight_array(x):
    """
    """
    x = x.ravel()
    x= x[x>0]

    ls = []
    for i in np.unique(x):
        n = len(x[x == i])/i
        ls.append([i]*int(n))

    ls = flatten(ls)    
    return np.array(ls)

def make_blob_veg(Nxcell, Nycell,sigma, fV):
    """
    """
    veg = np.random.uniform(0, 1, (Nxcell, Nycell)) 
    gauss = gaussian_filter(veg, sigma=sigma)
    threshold = np.percentile(gauss, fV*100)
    veg[gauss > threshold] = 0
    veg[gauss <= threshold] = 1
 
    return veg 

def autocorr(x):
    """
    autocorrelation using np
    """
    result = np.correlate(x, x, mode='full')
    return result[int(result.size/2):]/result.size


def fourier_autocorr(a):
    """
    autocorrelation using fourier
    """
    a = np.concatenate((a,np.zeros(len(a)-1))) # added zeros to your signal
    A = np.fft.fft(a)
    S = np.conj(A)*A
    c_fourier = np.fft.ifft(S)
    c_fourier = c_fourier[:(c_fourier.size//2)+1]
    return c_fourier/c_fourier.size/2



