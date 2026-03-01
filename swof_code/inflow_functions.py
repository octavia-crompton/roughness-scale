# coding=utf-8
""" gw_functions.py

Variable (Giraldez and Woolhiser's)
      L  (L) : hillslope length (m)
      alpha_eff : coefficient expressing surface conditions for the flow
      t_rain  : storm duration (s)
      rain  : rainfall intensity (m/s)
      Ks  :  saturated hydraulic conductivity (m/s)
      Ao  :   sorptivity, m/s^(1/2)
      t_pond  :  time of ponding (s)
"""
from __future__ import print_function

from read_SWOF import *
from plot_SWOF import *


import numpy as np
from scipy.optimize import fsolve

def compute_front(sim, thresh = 1e-5):
    """
    Returns the time at which the wetting front reaches each position on the slope

    Parameters:
    ----------
    sim : dict
        SVE simulation

    Returns
    ----------
    x_front : array_like
        x positions of the first characteristic (m)
    t_front : array_like
        t positions of the first characteristic (s)
    """
    front = []
    x_front = sim.yc.mean(0).copy()
    t_front = x_front.copy()

    for ind, x in enumerate(x_front):

        try:
            t_ind = np.where(sim['hc'].mean(1)[:, ind] > thresh)[0][0]
            t_r = sim.t[t_ind]
            t_front[ind] = t_r
        
        except IndexError:
            t_front[ind] = np.inf


    return x_front, t_front

def analytic_front(sim, characteristic = 1):
    """
    Compute position and time for wetting front, using kinematic assumption
    and continuity

    Parameters:
    ----------
    sim : pandas.Series
        SVE simulation

    characteristic : bool
        Use 1 if solution from method of characteristics is use

    Returns:
    -------
    d : list
        Front distance from divide (m)
    t : list
        Corresponding time (min)
    """
    x = sim.yc.mean(0)

    q_inflow = sim.bottom_imp_discharge/sim.L
    
    Kr = sim.So ** 0.5 / sim.alpha_v

    f = sim.Ks_v # Ksat (m/s)

    i = (sim.p - f)/3.6e5  #  m/s

    a = 2./3
    h_o = (q_inflow / Kr) ** (1./(a + 1))

    if characteristic:
        t = (1) / i  * (
        ( (i * x  + q_inflow )/ Kr) ** (1/(a+1)) - h_o)
    else: 
        t = (a + 1) / i  * (
            ((i * x  + q_inflow )/ Kr) ** (1/(a+1)) - h_o)
    return x, t 


################  Inflow plots ########################

def plot_bboundary(summary, out_dir,  dt = 30, ax=None, label = False):
    """
    """
    if not ax:
        fig, ax = plt.subplots(1, figsize = (7, 3))
    else:
        fig = plt.gcf()

    for key in summary.index:
        sim = summary.loc[key]
        
        b = read_boundary(sim, sim.key, out_dir)

        f = int(np.round(dt/np.unique(np.diff(b.t_boundary))[-1]))
        ax.plot(b.t_boundary[::f][1:], np.diff(b.bottom[::f]))

    ax.set_xlabel('minutes')
    ax.set_ylabel('m^3')

    return fig, ax


def plot_tboundary(summary,  out_dir, dt = 30, ax=None, label = False):
    """

    """
    if not ax:
        fig, ax = plt.subplots(1, figsize = (7, 3))
    else:
        fig = plt.gcf()

    for key in summary.index:
        sim = summary.loc[key]
        
        b = read_boundary(sim, sim.key, out_dir)

        f = int(np.round(dt/np.unique(np.diff(b.t_boundary))[-1]))
        ax.plot(b.t_boundary[::f][1:], np.diff(b.top[::f]))

    ax.set_xlabel('minutes')
    ax.set_ylabel(r'$m^3$')
    return fig, ax
    