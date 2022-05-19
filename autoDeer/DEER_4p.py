import importlib
import deerlab as dl
import numpy as np
import os,sys;
from scipy.integrate import cumulative_trapezoid

from autoDeer.File_Saving import save_file
import time
import importlib
import logging

log = logging.getLogger('core.DEER')


MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations


def IdentifyROI(dist:np.ndarray,r:np.ndarray,criterion:float = 0.99,plot:bool=False):
    """IdentifyROI _summary_

    Parameters
    ----------
    dist : np.ndarray
        The distance distribution.
    r : np.ndarray
        The distance axis
    criterion : float, optional
        What fraction of the distance distribution must be in the ROI, by default 0.99
    plot : bool, optional
        Plot the cumulative graphs, by default False
    """

    # Normlaize the distribution
    dist = dist / np.trapz(dist,r)
    # cumulative_dist = cumulative_trapezoid(dist,r,initial=0)
    # min_dist = r[np.argmin(np.abs(1 - cumulative_dist - criterion))]
    # max_dist = r[np.argmin(np.abs(cumulative_dist - criterion))]

    c_trapz_dist = np.zeros((dist.shape[0],dist.shape[0]))

    for i in range(0,dist.shape[0]):
        c_trapz_dist[i,i:] = cumulative_trapezoid(dist[i:],r[i:],initial=0)

    c_trapz_dist[(c_trapz_dist < criterion)] = 3
    ind = np.unravel_index(np.argmin(c_trapz_dist),c_trapz_dist.shape)
    min_dist = r[ind[0]]
    max_dist = r[ind[1]]

    # Enlarge ROI
    width = max_dist - min_dist
    max_dist = max_dist + width *0.25
    min_dist = min_dist - width *0.25

    if not plot:
        return[min_dist,max_dist]

    else:
        pass

def QualityControl(fit):

    if fit.regparam < 1:
        log.warning('Regularization Parameter is too small. Please try a different model')
        return 1
    else:
        return 0

    
