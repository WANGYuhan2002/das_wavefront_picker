# 
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import scipy.signal as signal
from scipy.ndimage import uniform_filter1d
import pandas as pd

#-- To import a function on a relative path:
import sys
sys.path.append("./")
from pydas_readers.readers import load_das_h5, write_das_h5
from pydas_readers.util import block_filters, block_cleaning

import os
import gc
from matplotlib import ticker
import matplotlib.dates as mdates


# !================= import the custom load function =================
import load_das_h5_PE2024 as load_das_h5_PE2024




# +
def sta_lta_centred(signal, nsta, nlta):
    """
    Calculates the ratio of the average of a^2 in a short-term (signal) window to a
    preceding long-term (noise) window. STA/LTA value is assigned to the end of the LTA
    /one sample before the start of the STA.

    Parameters
    ----------
    signal : array-like
        Signal array
    nsta : int
        Number of samples in short-term window
    nlta : int
        Number of samples in long-term window

    Returns
    -------
    sta / lta : array-like
        Ratio of a^2 in a short term average window to a preceding long term average
        window. STA/LTA value is assigned to end of LTA window / one sample before the
        start of STA window -- "centred".

    """
    # Cumulative sum to calculate moving average
    sta = np.cumsum(signal**2)
    sta = np.require(sta, dtype=float)
    lta = sta.copy()

    # Compute the STA and the LTA
    sta[nsta:] = sta[nsta:] - sta[:-nsta]
    sta[nsta:-nsta] = sta[nsta * 2 :]
    sta /= nsta

    lta[nlta:] = lta[nlta:] - lta[:-nlta]
    lta /= nlta

    sta[: (nlta - 1)] = 0
    sta[-nsta:] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(float).tiny
    idx = lta < dtiny
    lta[idx] = dtiny
    sta[idx] = 0.0

    # normalize the sta/lta by the max or MAD
    norm = 'MAD'
    if norm == 'MAD':
        mad = np.median(np.abs(sta/lta - np.median(sta/lta)))
        stalta_normed = sta / lta /mad
    if norm == 'MAX':
        stalta_normed = sta / lta /np.max(sta/lta) *2
    return stalta_normed

def sta_lta_centred_multichan(stream, nsta, nlta):
    """Apply STA/LTA to all channels in the stream."""
    csl_waterfall = np.zeros((stream.shape[0],stream.shape[1]))
    for idx_channel in range(0,stream.shape[1]):
        # calculate centered STA/LTA value, normalized by MAD, for each channel
        csl = sta_lta_centred(stream[:,idx_channel], nsta, nlta)
        csl_waterfall[:,idx_channel] = csl
    return csl_waterfall
# +
def trig_onechan_and_outindex(csl_onechan,threshold,minimal_event_interval,fs_das):
    """
    Trigger detection for a single channel.
    Returns indices of local maxima within triggered clusters.
    """
    csl = csl_onechan
    csl_01 = np.zeros(csl.shape)
    csl_01[csl > threshold] = 1
    index_1 = np.where(csl_01 == 1)
    
    # check if there is no event
    if len(index_1[0]) == 0:
        return np.array([])
    elif len(index_1[0]) > 0:
        # * Cluster the events with minimal_event_interval
        # index_2 saves the index of the first peak in each cluster
        index_2 = [] 
        index_2.append(index_1[0][0])
        for i in range(1, len(index_1[0])):
            if index_1[0][i] - index_1[0][i-1] > minimal_event_interval*fs_das:
                index_2.append(index_1[0][i])
        # print('start index',index_2)

        # index_3 saves the index of the last peak in each cluster
        index_3 = []
        index_3.append(index_1[0][-1])
        for i in range(len(index_1[0])-2, -1, -1):
            if index_1[0][i+1] - index_1[0][i] > minimal_event_interval*fs_das:
                index_3.append(index_1[0][i])
        index_3 = index_3[::-1]
        # print('end index' ,index_3)

        # pick the maximum in each cluster
        index_arg_max = []
        index_2 = np.array(index_2)
        index_3 = np.array(index_3)
        for i in range(len(index_2)):
            if index_2[i] == index_3[i]:
                index_arg_max.append(index_2[i])
            else:
                index_arg_max.append(index_2[i] + np.argmax(abs(csl[index_2[i]:index_3[i]])))
        return np.array(index_arg_max)
        
def trig_allchan_and_collapsed_index(csl,threshold,minimal_event_interval,fs_das):
    # // Legacy function: simple collapse without advanced tracking.
    # in index_t_trig, it is already sure that the (idx2-idx1)*fs_das>minimal_event_interval
    # then in index_arg_max, collapse the events that are too close to each other, and count how many are triggered
    index_arg_max = np.array([])
    score = np.array([])
    empty_no = 0
    csl_value = np.array([])
    
    for ichan in range(csl.shape[1]):
        index_t_trig = trig_onechan_and_outindex(csl[:,ichan],threshold,minimal_event_interval,fs_das)
        
        if len(index_t_trig) == 0:
            empty_no = empty_no + 1
            # print('------- empty in channel dist =',sel_dist[ichan])
        elif len(index_t_trig) > 0:
            if len(index_arg_max) == 0: # means still no initialization
                index_arg_max = index_t_trig
                #add a score to each index
                score = np.ones(len(index_t_trig))
                csl_value = csl[index_t_trig,ichan]
            else:
                for index_t in index_t_trig:
                    if np.min(abs(index_t - np.array(index_arg_max))) > minimal_event_interval*fs_das:
                        # print('new event','index_t',index_t)
                        index_arg_max = np.append(index_arg_max,index_t) # add a new index
                        score = np.append(score,1)
                        csl_value = np.append(csl_value,csl[index_t,ichan])
                    elif np.min(abs(index_t - np.array(index_arg_max))) < minimal_event_interval*fs_das:
                        idx_overlap = np.argmin(abs(index_t - np.array(index_arg_max)))
                        # print('overlap event','index_t',index_t)
                        # index_arg_max[idx_overlap] = index_t
                        score[idx_overlap] += 1
                        csl_value[idx_overlap] = csl_value[idx_overlap] + csl[index_t,ichan]
            # print(ichan,'index_arg_max',index_arg_max)
    if len(csl_value) == 0:
        return np.array([]), np.array([]), np.array([]), empty_no
    
    return index_arg_max, score, csl_value/score, empty_no
    
    
    
    


# !=============== v-1130 ========================start=========================


import numpy as np

def _merge_wavefronts(wavefronts, fs, dx, merge_params):
    """
    [Simplified] Smart Merging Function
    
    Logic:
    Judges purely based on the average time center of the wavefronts.
    If the average time difference between two wavefronts is less than 
    merge_params['merge_time_tol'], they are considered the same event and merged.
    
    Parameters
    ----------
    wavefronts : list
        List of dictionaries containing wavefront data.
    fs : float
        Sampling rate (Hz)
    dx : float
        Channel spacing (m) - (Not directly used in simplified version, kept for interface compatibility)
    merge_params : dict
        Must contain 'merge_time_tol' (seconds), defining the max average time diff for merging.
        
    Returns
    -------
    merged_wfs : list
        Merged list of wavefronts.
    """
    if not wavefronts:
        return []

    # 1. Get merge threshold (default 0.1s)
    time_tol_seconds = merge_params.get('merge_time_tol', 0.1) 
    time_tol_samples = time_tol_seconds * fs

    # 2. Sort by average time of wavefronts. This is crucial for the greedy merge algorithm.
    #    Ensures we always compare temporally adjacent wavefronts first.
    wavefronts.sort(key=lambda x: x['avg_time'])
    
    merged_wfs = []
    
    # Initialize the first wavefront as current candidate
    if len(wavefronts) > 0:
        current = wavefronts[0]
        
        for i in range(1, len(wavefronts)):
            candidate = wavefronts[i]
            
            # Calculate absolute difference in average time
            dt_avg = abs(candidate['avg_time'] - current['avg_time'])
            
            # --- Core Merge Logic ---
            if dt_avg < time_tol_samples:
                # Condition met, execute merge
                
                # 1. Merge details list (chan, time, val)
                combined_details = current['details'] + candidate['details']
                
                # 2. Deduplication: If the same channel exists in both, keep the earlier one.
                #    This aligns with the physical definition of "First Arrival".
                unique_dict = {}
                for d in combined_details:
                    ch = d[0]
                    t_idx = d[1]
                    if ch not in unique_dict:
                        unique_dict[ch] = d
                    else:
                        # Replace if new point arrives earlier
                        if t_idx < unique_dict[ch][1]:
                            unique_dict[ch] = d
                
                # Convert back to list and sort by channel ID for plotting/processing
                new_details = sorted(list(unique_dict.values()), key=lambda x: x[0])
                
                # 3. Extract new statistics
                w_times = np.array([p[1] for p in new_details])
                w_chans = np.array([p[0] for p in new_details])
                w_vals  = np.array([p[2] for p in new_details])
                
                # 4. Update current object properties
                current = {
                    'id': current['id'], # Keep ID constant or reset
                    'start_index': int(np.min(w_times)),
                    'end_index':   int(np.max(w_times)),
                    'duration_idx': int(np.max(w_times) - np.min(w_times)),
                    'avg_time':    np.mean(w_times), # Update avg_time, this affects next comparison (chain merge)
                    'score':       len(np.unique(w_chans)), # Update channel count score
                    'peak_val_sum': np.sum(w_vals),
                    'details':     new_details
                }
                # Note: We do not append 'candidate' to merged_wfs yet. We "swallow" it into 'current',
                # and continue comparing this updated 'current' with the next candidate.
                
            else:
                # Merge condition not met, 'current' is finalized.
                merged_wfs.append(current)
                # Set 'candidate' as the new 'current' for the next round.
                current = candidate
        
        # End of loop: Don't forget to append the last 'current'.
        merged_wfs.append(current)

    return merged_wfs


def trig_allchan_and_collapsed_index_and_save_wavefonts(
    csl, 
    threshold, 
    minimal_event_interval, 
    fs_das, 
    dx,                   
    tracking_params,      
    merge_params          
):
    """
    [V4.1 Fixed] Includes Slope/Acceleration Constraints + Simplified Merge Function
    
    Parameters
    ----------
    tracking_params : dict
        {
            'neighbor_tol': int,   # Spatial search range (channels)
            'time_tol': float,     # Time search window (s)
            'slope_tol': float     # Slope change tolerance (s/ch)
        }
    merge_params : dict
        {
            'merge_time_tol': float # [New] Allowed avg time diff for merging (s)
        }
    """
    n_time, n_chan = csl.shape
    
    # --- Unpack Parameters ---
    neighbor_tolerance = tracking_params.get('neighbor_tol', 3)
    neighbor_time_tolerance = tracking_params.get('time_tol', 0.05)
    
    # Get slope tolerance, None if not provided
    slope_tol_seconds = tracking_params.get('slope_tol', None)
    
    # Convert to sample units
    if slope_tol_seconds is not None:
        slope_tol_samples = slope_tol_seconds * fs_das 
    else:
        slope_tol_samples = None

    # --- A. Data Prep (Single Channel Triggering) ---
    triggers_by_chan = [[] for _ in range(n_chan)]
    empty_chan = 0
    
    # Iterate channels to find local triggers
    for ichan in range(n_chan):
        idxs = trig_onechan_and_outindex(csl[:, ichan], threshold, minimal_event_interval, fs_das)
        if len(idxs) == 0:
            empty_chan += 1
        else:
            vals = csl[idxs, ichan]
            for i, idx in enumerate(idxs):
                triggers_by_chan[ichan].append({
                    't_idx': idx, 'val': vals[i], 'used': False, 'ch_idx': ichan  
                })

    wavefronts = []
    
    # --- B. Chain Tracking (Core Logic) ---
    for start_ch in range(n_chan):
        for seed in triggers_by_chan[start_ch]:
            if seed['used']: continue
                
            current_chain = [seed]
            seed['used'] = True
            curr_node = seed
            
            # Initialize slope from previous step (Slope = dt / d_chan)
            prev_slope = None 
            
            while True:
                found_next = False
                
                # Collect all potential candidates in neighborhood
                all_candidates = []
                
                for gap in range(1, neighbor_tolerance + 1):
                    next_ch = curr_node['ch_idx'] + gap
                    if next_ch >= n_chan: break
                    
                    time_window_samples = int(neighbor_time_tolerance * fs_das * gap)
                    t_min = curr_node['t_idx'] - time_window_samples
                    t_max = curr_node['t_idx'] + time_window_samples
                    
                    for point in triggers_by_chan[next_ch]:
                        if point['used']: continue 
                        if t_min <= point['t_idx'] <= t_max:
                            # Cache candidate and calculate properties
                            dt = point['t_idx'] - curr_node['t_idx']
                            dch = gap # Spatial distance (in channels)
                            current_slope = dt / dch # Slope of this step
                            
                            all_candidates.append({
                                'point': point,
                                'dt_abs': abs(dt),
                                'slope': current_slope
                            })
                
                # --- Competition Mechanism ---
                best_next_wrapper = None
                
                if all_candidates:
                    # 1. Acceleration Constraint (Slope Consistency)
                    candidate_pool = all_candidates
                    if prev_slope is not None and slope_tol_samples is not None:
                        # Filter out points with sudden slope changes
                        candidate_pool = [
                            c for c in all_candidates 
                            if abs(c['slope'] - prev_slope) < slope_tol_samples
                        ]
                    
                    if candidate_pool:
                        # 2. Optimization Strategy
                        if prev_slope is not None:
                            # Choose the one with minimal slope change (Smoothest)
                            candidate_pool.sort(key=lambda c: abs(c['slope'] - prev_slope))
                        else:
                            # Choose the closest one in time
                            candidate_pool.sort(key=lambda c: c['dt_abs'])
                        
                        best_next_wrapper = candidate_pool[0]

                # --- Link Connection ---
                if best_next_wrapper:
                    best_next = best_next_wrapper['point']
                    best_next['used'] = True
                    current_chain.append(best_next)
                    
                    # Update slope memory
                    prev_slope = best_next_wrapper['slope']
                    
                    curr_node = best_next
                    found_next = True
                
                if not found_next: break
            
            # Encapsulate Wavefront
            # Only chains longer than 1 are considered valid wavefronts
            if len(current_chain) > 1:
                w_times = np.array([p['t_idx'] for p in current_chain])
                w_chans = np.array([p['ch_idx'] for p in current_chain])
                w_vals  = np.array([p['val'] for p in current_chain])
                
                wf = {
                    'id': len(wavefronts) + 1,
                    'start_index': int(np.min(w_times)),
                    'end_index':   int(np.max(w_times)),
                    'duration_idx':int(np.max(w_times) - np.min(w_times)),
                    'avg_time':    np.mean(w_times),
                    'score':       len(np.unique(w_chans)),
                    'peak_val_sum':np.sum(w_vals),
                    'details': sorted(zip(w_chans.astype(int), w_times.astype(int), w_vals), key=lambda x: x[0])
                }
                wavefronts.append(wf)

    # --- C. Call Simplified Merge Function ---
    if len(wavefronts) > 0:
        wavefronts = _merge_wavefronts(wavefronts, fs_das, dx, merge_params)

    # --- D. Return Results ---
    if not wavefronts:
        return np.array([]), np.array([]), np.array([]), np.array([]), empty_chan, []

    # Final Sort to ensure output order
    wavefronts.sort(key=lambda x: x['start_index'])
    
    index_start = np.array([wf['start_index'] for wf in wavefronts])
    index_end   = np.array([wf['end_index'] for wf in wavefronts])
    scores      = np.array([wf['score'] for wf in wavefronts])
    avg_vals    = np.array([wf['peak_val_sum']/wf['score'] for wf in wavefronts])
    
    return index_start, index_end, scores, avg_vals, empty_chan, wavefronts

#!=============== v-1130 ========================end=========================