#! /usr/bin/env python

# ============================ import libraries ============================
# built-in libraries
import os
import csv
import math
import logging
import pathlib
import platform
import itertools
from functools import partial
from collections import Counter, defaultdict
from multiprocessing import Queue, Process, Pool
from copy import copy

# third-party libraries
import numpy as np
import pandas as pd
from pandas.io import pickle
from pyteomics import mzml
import opentims_bruker_bridge # pip install opentims_bruker_bridge
from opentimspy.opentims import OpenTIMS # pip install opentimspy


# ============================ define functions ============================
# data preparation
def normalize_path(path):
    return os.path.normpath(os.path.abspath(path))

def prepare(msfile, library, outpath):
    msfile = normalize_path(msfile)
    if not os.path.exists(msfile):
        raise FileNotFoundError(msfile)
    if not (msfile.endswith(".d") or msfile.endswith(".raw") or msfile.endswith(".mzML")):
        raise ValueError("Only MS file in '.d', '.raw' or '.mzML' format are supported!")
    library = normalize_path(library)
    if not os.path.exists(library):
        raise FileNotFoundError(library)
    outpath = normalize_path(outpath)
    return msfile, library, outpath

# convert raw into mzml
def convert_raw_to_mzml(rawfile):
    if platform.system() == "Windows":
        RawFileParser = normalize_path(f"{os.getcwd()}/third_party/ThermoRawFileParser.exe")
        if not os.path.exists(RawFileParser):
            raise FileNotFoundError(RawFileParser)
        convert_status = os.system("%s -i=%s -o=%s -f=1 > NUL 2> NUL" % (RawFileParser, rawfile, os.path.dirname(rawfile)))
    elif platform.system() == "Linux":
        RawFileParser = normalize_path(f"{os.getcwd()}/third_party/ThermoRawFileParser")
        if not os.path.exists(RawFileParser):
            raise FileNotFoundError(RawFileParser)
        convert_status = os.system("%s -i=%s -o=%s -f=1 > /dev/null 2>&1" % (RawFileParser, rawfile, os.path.dirname(rawfile)))
    else:
        raise Exception(f"Unsupported platform: {platform.system()}")
    return convert_status
    
# extract MS1 spectra
def extract_ms1(msfile):
    if msfile.endswith(".d"):
        path = pathlib.Path(msfile)
        D = OpenTIMS(path)
        ms1_frames = D.query_iter(D.ms1_frames, columns = ('frame', 'retention_time', 'mz', 'intensity', 'inv_ion_mobility'))
        ms1 = []
        for idx, frame in enumerate(ms1_frames):
            new_frame = {'scan idx':idx + 1, 
                        'retention time':round(frame['retention_time'][0]/60, 4), 
                        'm/z array':frame['mz'], 
                        'intensity array':frame['intensity'],
                        'mean inverse reduced ion mobility array':frame['inv_ion_mobility']}   
            ms1.append(new_frame)
    elif msfile.endswith(".mzML"):
        spectra = mzml.read(msfile)
        ms1 = []
        idx = 0
        for spectrum in spectra:
            if spectrum['ms level'] == 1:
                idx += 1
                new_spectrum = {'scan idx':idx, 
                                'retention time':round(spectrum['scanList']['scan'][0]['scan start time'], 4), 
                                'm/z array':spectrum['m/z array'], 
                                'intensity array':spectrum['intensity array']}
                if 'mean inverse reduced ion mobility array' in spectrum:
                    new_spectrum['mean inverse reduced ion mobility array'] = spectrum['mean inverse reduced ion mobility array']
                ms1.append(new_spectrum)
    return ms1

# process MS1 spectra
def combine_ms1(ms1):
    """
    func: combine peaks with the same m/z in each ms1 spectrum
    input: ms1 (with ion mobility)
    output: ms1 (without ion mobility)
    """
    new_ms1 = []
    for z in ms1:
        mz_array = z['m/z array']
        intensity_array = z['intensity array']
        mz_intensity = {}
        for i in range(len(mz_array)):
            if mz_array[i] in mz_intensity:
                mz_intensity[mz_array[i]] += intensity_array[i]
            else:
                mz_intensity[mz_array[i]] = intensity_array[i]
        sorted_mz_intensity = dict(sorted(mz_intensity.items(), key = lambda item: item[0]))
        new_mz_array = np.array([mz for mz in sorted_mz_intensity.keys()])
        new_intensity_array = np.array([intv for intv in sorted_mz_intensity.values()])
        z['m/z array'] = new_mz_array
        z['intensity array'] = new_intensity_array
        del z['mean inverse reduced ion mobility array']
        new_ms1.append(z)
    return new_ms1

def filter_ms1(ms1, min_mz, max_mz, min_intensity):
    """
    func: filter peaks in each spectrum out of range
    """
    new_ms1 = []
    for z in ms1:
        mz_array = z['m/z array']
        intensity_array = z['intensity array']
        peaks_range =  (mz_array >= min_mz) & (mz_array <= max_mz) & (intensity_array >= min_intensity)
        z['m/z array'] = mz_array[peaks_range]
        z['intensity array'] = intensity_array[peaks_range]
        new_ms1.append(z)
    return new_ms1

def centroid_ms1(ms1, mass_diff):
    """
    func: convert ms1 spectra into centroid mode
    """
    new_ms1 = []
    for z in ms1:
        best_mz, best_int = 0, 0
        prev_mz, prev_int = False, False
        mzs, ints = [], []
        for cur_mz, cur_int in zip(z['m/z array'], z['intensity array']):
            if prev_mz is False:
                best_mz = cur_mz
                best_int = cur_int
            elif cur_mz - prev_mz > mass_diff:
                mzs.append(best_mz)
                ints.append(best_int)
                best_mz = cur_mz
                best_int = cur_int
            elif cur_mz - prev_mz <= mass_diff:
                if best_int > prev_int and cur_int > prev_int:
                    mzs.append(best_mz)
                    ints.append(best_int)
                    best_mz = cur_mz
                    best_int = cur_int
                elif cur_int > best_int:
                    best_mz = cur_mz
                    best_int = cur_int
            prev_mz = cur_mz
            prev_int = cur_int
        mzs.append(best_mz)
        ints.append(best_int)
        z['m/z array'] = np.array(mzs)
        z['intensity array'] = np.array(ints)
        new_ms1.append(z)
    return new_ms1

def process_ms1(ms1, min_mz, max_mz, min_intensity, profile, mass_diff):
    if all('mean inverse reduced ion mobility array' in z for z in ms1):
        ms1 = combine_ms1(ms1)
    ms1 = filter_ms1(ms1, min_mz, max_mz, min_intensity)
    max_mz_value = max((z['m/z array'].max() for z in ms1), default = 0)
    if profile:
        ms1 = centroid_ms1(ms1, mass_diff)
    return ms1, max_mz_value

# build hills
def get_fast_array(array, step, idx):
    fast_dict = {}
    fast_array = list((array / step).astype(int))
    for idx, fm in zip(idx, fast_array):
        if fm not in fast_dict:
            fast_dict[fm] = [idx, ]
        else:
            fast_dict[fm].append(idx)
    return fast_array, fast_dict

def build_hills(ms1, hill_mass_tol, mz_step, im_tol):
    hills_dict = {'hill_idx_array':[], 
                  'orig_idx_array':[], 
                  'scan_idx_array':[], 
                  'mz_array':[], 
                  'intensity_array':[]}
    prev_fast_mz_dict = {}
    if im_tol > 0:
        hills_dict['im_array'] = []
        prev_fast_im_dict = {}
    last_idx = -1
    prev_idx = -1
    for z_idx, z in enumerate(ms1):
        len_mz_array = len(z['m/z array'])
        hills_dict['hill_idx_array'].extend(list(range(last_idx+1, last_idx+1+len_mz_array, 1)))
        hills_dict['orig_idx_array'].extend(range(len_mz_array))
        hills_dict['scan_idx_array'].extend([z_idx] * len_mz_array)
        hills_dict['mz_array'].extend(z['m/z array'])
        hills_dict['intensity_array'].extend(z['intensity array'])
        sorted_idx = np.argsort(z['intensity array'])[::-1]
        sorted_mz_array = z['m/z array'][sorted_idx]
        fast_mz_array, fast_mz_dict = get_fast_array(sorted_mz_array, mz_step, sorted_idx)
        if im_tol > 0:
            hills_dict['im_array'].extend(z['mean inverse reduced ion mobility array'])
            sorted_im_array = z['mean inverse reduced ion mobility array'][sorted_idx]
            fast_im_array, fast_im_dict = get_fast_array(sorted_im_array, im_tol, sorted_idx)
        banned_prev_idx_set = set()
        for idx, fm, fi in zip(list(sorted_idx), fast_mz_array, (fast_mz_array if im_tol == 0 else fast_im_array)):
            flag1 = fm in prev_fast_mz_dict
            flag2 = fm-1 in prev_fast_mz_dict
            flag3 = fm+1 in prev_fast_mz_dict
            if flag1 or flag2 or flag3:
                if flag1:
                    all_idx = prev_fast_mz_dict[fm]
                    if flag2:
                        all_idx += prev_fast_mz_dict[fm-1]
                    if flag3:
                        all_idx += prev_fast_mz_dict[fm+1]
                elif flag2:
                    all_idx = prev_fast_mz_dict[fm-1]
                    if flag3:
                        all_idx += prev_fast_mz_dict[fm+1]
                elif flag3:
                    all_idx = prev_fast_mz_dict[fm+1]
                if im_tol > 0:
                    flag1_im = fi in prev_fast_im_dict
                    flag2_im = fi-1 in prev_fast_im_dict
                    flag3_im = fi+1 in prev_fast_im_dict
                best_intensity = 0
                best_idx_prev = 0
                mz_cur = z['m/z array'][idx]
                all_prevs = [[idx_prev, ms1[z_idx-1]['intensity array'][idx_prev]] for idx_prev in all_idx if (idx_prev not in banned_prev_idx_set and (im_tol == 0 or (flag2_im and idx_prev in prev_fast_im_dict[fi-1]) or (flag1_im and idx_prev in prev_fast_im_dict[fi]) or (flag3_im and idx_prev in prev_fast_im_dict[fi+1])))]
                for idx_prev, cur_intensity in all_prevs:
                    cur_mass_diff = abs((mz_cur - ms1[z_idx-1]['m/z array'][idx_prev]) / mz_cur * 1e6)
                    if cur_mass_diff <= hill_mass_tol and cur_intensity >= best_intensity:
                        best_intensity = cur_intensity
                        best_idx_prev = idx_prev
                        hills_dict['hill_idx_array'][last_idx+1+idx] = hills_dict['hill_idx_array'][prev_idx+1+idx_prev]
                if best_idx_prev != 0:
                    banned_prev_idx_set.add(best_idx_prev)
        prev_fast_mz_dict = fast_mz_dict
        if im_tol > 0:
            prev_fast_im_dict = fast_im_dict
        prev_idx = last_idx
        last_idx = last_idx + len_mz_array
    return hills_dict

# split hills
def split_peaks_python(qout, hills_dict, data_for_analyse_tmp, counter_hills_idx, sorted_idx_child_process, sorted_idx_array_child_process, i, checked_id, min_length_hill, hillValleyFactor, win_sys = False):
    new_index_list = split_peaks(hills_dict, data_for_analyse_tmp, counter_hills_idx, sorted_idx_child_process, sorted_idx_array_child_process, i, checked_id, min_length_hill, hillValleyFactor)
    if win_sys:
        return (i, list(new_index_list))
    else:
        qout.put((i, list(new_index_list)))
        qout.put(None)

def split_peaks(hills_dict, data_for_analyse_tmp, counter_hills_idx, sorted_idx_child_process, sorted_idx_array_child_process, nproc, checked_id, min_length_hill, hillValleyFactor):
    new_index_list = copy(sorted_idx_array_child_process)
    idx_start = 0
    idx_end = 0
    cur_new_idx = checked_id + max(sorted_idx_child_process) + 1
    for hill_idx in sorted_idx_child_process:
        hill_length = counter_hills_idx[hill_idx]
        idx_end = idx_start + hill_length
        if hill_length >= min_length_hill * 2:
            tmp_scans = hills_dict['scan_idx_array'][checked_id+idx_start:checked_id+idx_end]
            tmp_orig_idx = hills_dict['orig_idx_array'][checked_id+idx_start:checked_id+idx_end]
            tmp_intensity = [data_for_analyse_tmp[scan_val][orig_idx_val] for orig_idx_val, scan_val in zip(tmp_orig_idx, tmp_scans)]
            smothed_intensity = meanfilt(tmp_intensity, 3)
            c_len = hill_length - min_length_hill
            idx = int(min_length_hill) - 1
            min_idx_list = []
            min_val = 0
            l_idx = 0
            recheck_r_r = []
            while idx <= c_len:
                if len(min_idx_list) and idx >= min_idx_list[-1] + min_length_hill:
                    l_idx = min_idx_list[-1]
                l_r = max(smothed_intensity[l_idx:idx]) / float(smothed_intensity[idx])
                if l_r >= hillValleyFactor:
                    r_r = max(smothed_intensity[idx+1:]) / float(smothed_intensity[idx])
                    if r_r >= hillValleyFactor:
                        mult_val = l_r * r_r
                        include_factor = (1 if l_r > r_r else 0)
                        if (min_length_hill <= idx + include_factor <= c_len):
                            if not len(min_idx_list) or idx + include_factor >= min_idx_list[-1] + min_length_hill:
                                min_idx_list.append(idx + include_factor)
                                recheck_r_r.append(idx)
                                min_val = mult_val
                            elif mult_val > min_val:
                                min_idx_list[-1] = idx + include_factor
                                recheck_r_r[-1] = idx
                                min_val = mult_val
                idx += 1
            if len(min_idx_list):
                for min_idx, end_idx, recheck_idx in zip(min_idx_list, min_idx_list[1:] + [idx_start+hill_length, ], recheck_r_r):
                    r_r = max(smothed_intensity[recheck_idx+1:end_idx]) / float(smothed_intensity[recheck_idx])
                    if r_r >= hillValleyFactor:
                        new_index_list[idx_start+min_idx:idx_start+hill_length] = cur_new_idx
                        cur_new_idx += 1
        idx_start = idx_end
    return new_index_list
    
def meanfilt(data, window_width):
    cumsum_vec = np.cumsum(data, dtype=float)
    cumsum_vec[window_width:] = cumsum_vec[window_width:] - cumsum_vec[:-window_width]
    ma_vec = data[:1] + list(cumsum_vec / window_width) + data[-1:]
    return ma_vec

def split_peaks_multi(hills_dict, spectra, hill_min_length, hillValleyFactor, n_procs):
    hills_dict['orig_idx_array'] = np.array(hills_dict['orig_idx_array'])
    hills_dict['scan_idx_array'] = np.array(hills_dict['scan_idx_array'])
    hills_dict['hill_idx_array'] = np.array(hills_dict['hill_idx_array'])
    counter_hills_idx = dict(Counter(hills_dict['hill_idx_array']))
    tmp_hill_length = np.array([counter_hills_idx[hill_idx] for hill_idx in hills_dict['hill_idx_array']])
    idx_range = tmp_hill_length >= hill_min_length
    hills_dict['hill_idx_array'] = hills_dict['hill_idx_array'][idx_range]
    hills_dict['scan_idx_array'] = hills_dict['scan_idx_array'][idx_range]
    hills_dict['orig_idx_array'] = hills_dict['orig_idx_array'][idx_range]
    all_sets = []
    all_sorted_idx = []
    if len(hills_dict['orig_idx_array']):
        idx_sort = np.argsort(hills_dict['hill_idx_array'] + ((hills_dict['scan_idx_array'] + 1) / (hills_dict['scan_idx_array'].max() + 2)))
        hills_dict['hill_idx_array'] = hills_dict['hill_idx_array'][idx_sort]
        hills_dict['scan_idx_array'] = hills_dict['scan_idx_array'][idx_sort]
        hills_dict['orig_idx_array'] = hills_dict['orig_idx_array'][idx_sort]
        hills_dict['hill_idx_array_unique'] = sorted(list(set(hills_dict['hill_idx_array'])))
        intensities_list = [spectrum['intensity array'] for spectrum in spectra]
        if n_procs == 1:
            qout = []
            procs = []
            new_idx_res = dict()
            checked_id = 0
            sorted_idx_child_process = list(hills_dict['hill_idx_array_unique'])
            idx_unique_set = set(sorted_idx_child_process)
            local_idx = np.array([z in idx_unique_set for z in list(hills_dict['hill_idx_array'])])
            sorted_idx_array_child_process = hills_dict['hill_idx_array'][local_idx]
            sorted_idx_child_process = sorted(list(idx_unique_set))
            i = 0
            qout = split_peaks_python(qout, hills_dict, intensities_list, counter_hills_idx, sorted_idx_child_process, sorted_idx_array_child_process, i, checked_id, hill_min_length, hillValleyFactor, win_sys = True)
            new_idx_res[qout[0]] = qout[1]
        else:
            qout = Queue()
            new_idx_res = dict()
            procs = []
            ar2 = []
            len_full = len(hills_dict['hill_idx_array_unique'])
            if len_full <= 1000 * n_procs:
                n_procs = 1
            step = int(len_full / n_procs) + 1
            checked_id = 0
            for i in range(n_procs):
                sorted_idx_child_process = list(hills_dict['hill_idx_array_unique'][i*step:i*step+step])
                idx_unique_set = set(sorted_idx_child_process)
                all_sets.append(idx_unique_set)
                local_idx = np.array([z in idx_unique_set for z in list(hills_dict['hill_idx_array'])])
                sorted_idx_array_child_process = hills_dict['hill_idx_array'][local_idx]
                ar2.append(sorted_idx_array_child_process)
                sorted_idx_child_process = sorted(list(idx_unique_set))
                all_sorted_idx.append(local_idx)
                p = Process(
                    target=split_peaks_python,
                    args=(qout, hills_dict, intensities_list, counter_hills_idx, sorted_idx_child_process, sorted_idx_array_child_process, i, checked_id, hill_min_length, hillValleyFactor))
                checked_id += len(sorted_idx_array_child_process)
                p.start()
                procs.append(p)
            for _ in range(n_procs):
                for ready_child_process in iter(qout.get, None):
                    new_idx_res[ready_child_process[0]] = ready_child_process[1]
            for p in procs:
                p.join()
    final_idx_array = []
    last_id = 1
    for i in range(n_procs):
        added_idx_map = {}
        for idx_val in new_idx_res[i]:
            if idx_val not in added_idx_map:
                added_idx_map[idx_val] = int(last_id)
                last_id += 1
            final_idx_array.append(added_idx_map[idx_val])
    hills_dict['hill_idx_array'] = list(final_idx_array)
    del hills_dict['hill_idx_array_unique']
    return hills_dict

# process hills
def process_hills(hills_dict, data_for_analyse_tmp, mz_step, paseftol, dia = False):
    counter_hills_idx = Counter(hills_dict['hill_idx_array'])
    if dia is False:
        min_length_hill = 2
    else:
        min_length_hill = 1
    hills_dict['orig_idx_array'] = np.array(hills_dict['orig_idx_array'])
    hills_dict['scan_idx_array'] = np.array(hills_dict['scan_idx_array'])
    hills_dict['hill_idx_array'] = np.array(hills_dict['hill_idx_array'])
    tmp_hill_length = np.array([counter_hills_idx[hill_idx] for hill_idx in hills_dict['hill_idx_array']])
    idx_minl = tmp_hill_length >= min_length_hill
    hills_dict['hill_idx_array'] = hills_dict['hill_idx_array'][idx_minl]
    hills_dict['scan_idx_array'] = hills_dict['scan_idx_array'][idx_minl]
    hills_dict['orig_idx_array'] = hills_dict['orig_idx_array'][idx_minl]
    if len(hills_dict['hill_idx_array']):
        idx_sort = np.argsort(hills_dict['hill_idx_array'] + (hills_dict['scan_idx_array'] / (hills_dict['scan_idx_array'].max()+1)))
        hills_dict['hill_idx_array'] = hills_dict['hill_idx_array'][idx_sort]
        hills_dict['scan_idx_array'] = hills_dict['scan_idx_array'][idx_sort]
        hills_dict['orig_idx_array'] = hills_dict['orig_idx_array'][idx_sort]
        hills_dict['hill_idx_array_unique'] = sorted(list(set(hills_dict['hill_idx_array'])))
        hills_dict['hill_mz_median'] = []
        hills_dict['hill_mz_median_fast_dict'] = defaultdict(list)
        if paseftol > 0:
            hills_dict['hill_im_median'] = []
            hills_dict['hill_im_median_fast_dict'] = defaultdict(set)
        hills_dict['hill_intensity_array'] = []
        hills_dict['hill_scan_set'] = []
        hills_dict['hill_scan_list'] = []
        hills_dict['hill_lengths'] = []
        hills_dict['tmp_mz_array'] = []
        hills_dict['scan_idx_array'] = list(hills_dict['scan_idx_array'])
        hills_dict['orig_idx_array'] = list(hills_dict['orig_idx_array'])
        hills_dict['hill_idx_array'] = list(hills_dict['hill_idx_array'])
        idx_start = 0
        idx_end = 0
        for idx_1, hill_idx in enumerate(hills_dict['hill_idx_array_unique']):
            hill_length = counter_hills_idx[hill_idx]
            idx_end = idx_start + hill_length
            tmp_scans = hills_dict['scan_idx_array'][idx_start:idx_end]
            tmp_orig_idx = hills_dict['orig_idx_array'][idx_start:idx_end]
            tmp_intensity = [data_for_analyse_tmp[scan_val]['intensity array'][orig_idx_val] for orig_idx_val, scan_val in zip(tmp_orig_idx, tmp_scans)]
            tmp_mz_array = [data_for_analyse_tmp[scan_val]['m/z array'][orig_idx_val] for orig_idx_val, scan_val in zip(tmp_orig_idx, tmp_scans)]
            mz_median = 0
            i_sum_tmp = 0
            for mz_val_tmp, i_val_tmp in zip(tmp_mz_array, tmp_intensity):
                mz_median += mz_val_tmp * i_val_tmp
                i_sum_tmp += i_val_tmp
            mz_median = mz_median / i_sum_tmp
            if paseftol > 0:
                tmp_im_array = [data_for_analyse_tmp[scan_val]['mean inverse reduced ion mobility array'][orig_idx_val] for orig_idx_val, scan_val in zip(tmp_orig_idx, tmp_scans)]
                im_median = np.average(tmp_im_array, weights=tmp_intensity)
            tmp_scans_list = tmp_scans
            tmp_scans_set = set(tmp_scans)
            idx_start = idx_end
            hills_dict['hill_mz_median'].append(mz_median)
            mz_median_int = int(mz_median/mz_step)
            tmp_val = (idx_1, tmp_scans_list[0], tmp_scans_list[-1])
            hills_dict['hill_mz_median_fast_dict'][mz_median_int-1].append(tmp_val)
            hills_dict['hill_mz_median_fast_dict'][mz_median_int].append(tmp_val)
            hills_dict['hill_mz_median_fast_dict'][mz_median_int+1].append(tmp_val)
            if paseftol > 0:
                hills_dict['hill_im_median'].append(im_median)
                im_median_int = int(im_median/paseftol)
                hills_dict['hill_im_median_fast_dict'][im_median_int-1].add(idx_1)
                hills_dict['hill_im_median_fast_dict'][im_median_int].add(idx_1)
                hills_dict['hill_im_median_fast_dict'][im_median_int+1].add(idx_1)
            hills_dict['hill_intensity_array'].append(tmp_intensity)
            hills_dict['hill_scan_set'].append(tmp_scans_set)
            hills_dict['hill_scan_list'].append(tmp_scans_list)
            hills_dict['hill_lengths'].append(hill_length)
            hills_dict['tmp_mz_array'].append(tmp_mz_array)
        hills_dict['hill_intensity_apex'] = [None] * len(hills_dict['hill_idx_array_unique'])
        hills_dict['hill_scan_apex'] = [None] * len(hills_dict['hill_idx_array_unique'])
    return hills_dict

def get_and_calc_apex_intensity_and_scan(hills_dict, idx_1):
    hill_intensity_apex_1 = hills_dict['hill_intensity_apex'][idx_1]
    hill_scan_apex_1 = hills_dict['hill_scan_apex'][idx_1]
    if hill_intensity_apex_1 is None:
        hill_intensity_apex_1 = 0
        for int_val, scan_val in zip(hills_dict['hill_intensity_array'][idx_1], hills_dict['hill_scan_list'][idx_1]):
            if int_val > hill_intensity_apex_1:
                hill_intensity_apex_1 = int_val
                hill_scan_apex_1 = scan_val
        hills_dict['hill_intensity_apex'][idx_1] = hill_intensity_apex_1
        hills_dict['hill_scan_apex'][idx_1] = hill_scan_apex_1
    return hills_dict, hill_intensity_apex_1, hill_scan_apex_1

def process_hills_extra(hills_dict, RT_dict):
    hills_features = {}
    for idx_1 in range(len(hills_dict['hill_idx_array_unique'])):
        hill_feature = {}
        hills_dict, hill_intensity_apex_1, hill_scan_apex_1 = get_and_calc_apex_intensity_and_scan(hills_dict, idx_1)
        hill_feature['mz'] = hills_dict['hill_mz_median'][idx_1]
        hill_feature['nScans'] = hills_dict['hill_lengths'][idx_1]
        hill_feature['rtApex'] = RT_dict[hill_scan_apex_1]
        hill_feature['scanApex'] = hill_scan_apex_1
        hill_feature['intensityApex'] = hill_intensity_apex_1
        hill_feature['intensitySum'] = sum(hills_dict['hill_intensity_array'][idx_1])
        hill_feature['rtStart'] = RT_dict[hills_dict['hill_scan_list'][idx_1][0]]
        hill_feature['rtEnd'] = RT_dict[hills_dict['hill_scan_list'][idx_1][-1]]
        if 'hill_im_median' in hills_dict:
            hill_feature['im'] = hills_dict['hill_im_median'][idx_1]
        else:
            hill_feature['im'] = 0
        hill_feature['hill_scan_list'] = hills_dict['hill_scan_list'][idx_1]
        hill_feature['hill_rt_list'] = [RT_dict[scan] for scan in hills_dict['hill_scan_list'][idx_1]]
        hill_feature['hill_intensity_list'] = hills_dict['hill_intensity_array'][idx_1]
        hills_features[idx_1] = hill_feature
    return hills_dict, hills_features

# detect BrIS envelopes
def get_candidate_hills_features(hills_features, mz, mz_tol):
    mz_range = [mz - mz * mz_tol * 1e-6, mz + mz * mz_tol * 1e-6]
    candidate_hills_features = {idx:features for idx, features in hills_features.items() if mz_range[0] <= features['mz'] <= mz_range[1]}
    if len(candidate_hills_features):
        candidate_hills_features = dict(sorted(candidate_hills_features.items(), key = lambda item:item[1]['intensityApex'], reverse = True))
    return candidate_hills_features

def cal_dot_product(array1, array2):
    return sum(i*j for i, j in zip(array1, array2))

def cal_mold_length(array):
    return math.sqrt(sum(i**2 for i in array))

def detect_BrIS_envolopes(library_dict, proton_mass, hills_features, mass_tolerance, nIsotopes, min_correlation, nFeatures):
    banned_hills_idx = []
    all_BrISs_envolopes = {}
    for BrIS_accession, info in library_dict.items():
        isotopes_mass = [float(x) for x in info['isotopes mass'].split(';')]
        isotopes_intensity_distribution = [float(x) for x in info['isotopes intensity distribution'].split(';')]
        BrIS_envolopes = []
        for charge in [2, 3]:
            iso0_mz = (isotopes_mass[0] + charge * proton_mass) / charge
            iso0_candidate_hills_features = get_candidate_hills_features(hills_features, iso0_mz, mass_tolerance)
            if len(iso0_candidate_hills_features):
                tmp_charge_envolopes = []
                for iso0_hill_idx, iso0_hill_features in iso0_candidate_hills_features.items():
                    iso0_hill_mz = iso0_hill_features['mz']
                    iso0_hill_scans = iso0_hill_features['hill_scan_list']
                    iso0_hill_intensities = iso0_hill_features['hill_intensity_list']
                    nIsotopes_candidates = []
                    for i in range(1, nIsotopes):
                        isoi_mz = iso0_hill_mz + (isotopes_mass[i] - isotopes_mass[0]) / charge
                        isoi_candidate_hills_features = get_candidate_hills_features(hills_features, isoi_mz, mass_tolerance)
                        isoi_candidate_hills_features = {idx:features for idx, features in isoi_candidate_hills_features.items() if idx not in banned_hills_idx}
                        if len(isoi_candidate_hills_features):
                            isoi_candidates = []
                            for isoi_hill_idx, isoi_hill_features in isoi_candidate_hills_features.items():
                                isoi_hill_scans = isoi_hill_features['hill_scan_list']
                                isoi_hill_intensities = isoi_hill_features['hill_intensity_list']
                                scan_intersection = list(set(iso0_hill_scans) & set(isoi_hill_scans))
                                if len(scan_intersection):
                                    iso0_intensities = [intv for scan, intv in zip(iso0_hill_scans, iso0_hill_intensities) if scan in scan_intersection]
                                    isoi_intensities = [intv for scan, intv in zip(isoi_hill_scans, isoi_hill_intensities) if scan in scan_intersection]
                                    dot_product = cal_dot_product(np.array(iso0_intensities).astype(np.float64), np.array(isoi_intensities).astype(np.float64))
                                    iso0_hill_intensities_mold_length = cal_mold_length(np.array(iso0_hill_intensities).astype(np.float64))
                                    isoi_hill_intensities_mold_length = cal_mold_length(np.array(isoi_hill_intensities).astype(np.float64))
                                    rt_cos_cor = dot_product / (iso0_hill_intensities_mold_length * isoi_hill_intensities_mold_length)
                                    if rt_cos_cor >= min_correlation:
                                        tmp_isotope_dict = {'isotope_number':i, 
                                                            'isotope_hill_idx':isoi_hill_idx, 
                                                            'isotope_hill_features':isoi_hill_features, 
                                                            'rt_cos_cor':rt_cos_cor}
                                        isoi_candidates.append(tmp_isotope_dict)
                            if len(isoi_candidates):
                                nIsotopes_candidates.append(isoi_candidates)
                        if len(nIsotopes_candidates) < i:
                            break
                    if len(nIsotopes_candidates) == nIsotopes - 1:
                        iso0_hill_intensityApex = iso0_hill_features['intensityApex']
                        all_isos_theoretical_intensityApex = [(isotopes_intensity_distribution[i] / isotopes_intensity_distribution[0]) * iso0_hill_intensityApex for i in range(nIsotopes)]
                        best_cor = 0
                        best_com = []
                        for iter_candidates in itertools.product(*nIsotopes_candidates):
                            all_isos_experimental_intensityApex = [iso0_hill_intensityApex] + [iso['isotope_hill_features']['intensityApex'] for iso in iter_candidates]
                            intensityApex_cos_cor = np.corrcoef(all_isos_theoretical_intensityApex, all_isos_experimental_intensityApex)[0,1]
                            if intensityApex_cos_cor >= max(min_correlation, best_cor):
                                best_cor = intensityApex_cos_cor
                                best_com = list(iter_candidates)
                                intensityApex_for_cos_cor = [all_isos_theoretical_intensityApex, all_isos_experimental_intensityApex]
                        if best_cor:
                            isos_hill_idx = [iso0_hill_idx] + [iso['isotope_hill_idx'] for iso in best_com]
                            tmp_all_isotopes_dict = {'iso0':iso0_hill_features,}
                            for i, iso in enumerate(best_com):
                                tmp_all_isotopes_dict["iso%d" % (i+1)] = iso['isotope_hill_features']
                            iso0_rtApex = iso0_hill_features['rtApex']
                            tmp_all_isotopes_dict.update({'iso0_rtApex':iso0_rtApex,
                                                          'iso0_intensityApex':iso0_hill_intensityApex,
                                                          'charge':charge,
                                                          'isotopes':best_com,
                                                          'intensityApex_cos_cor':best_cor,
                                                          'intensityApex_for_cos_cor':intensityApex_for_cos_cor})
                            tmp_charge_envolopes.append(tmp_all_isotopes_dict)
                            banned_hills_idx.extend(isos_hill_idx)      
                    if len(tmp_charge_envolopes) >= nFeatures:
                        break
                if len(tmp_charge_envolopes):
                    BrIS_envolopes.extend(tmp_charge_envolopes)
        if len(BrIS_envolopes):
            all_BrISs_envolopes[BrIS_accession] = BrIS_envolopes
    return all_BrISs_envolopes

# verify BrIS 
def calculate_sublists(list1):
    all_sublists = []
    for i in range(len(list1) + 1):
        sublists = itertools.combinations(list1, i)
        all_sublists.extend(sublists)
    return all_sublists
    
def process_candidates(candidates, BrIS_ready, BrIS_ready_refRT, refRT_ordered_BrIS_ready_sublists):
    canRT, intensity = candidates
    BrIS_ready_canRT = dict(zip(BrIS_ready, canRT))
    BrIS_ready_intensity = dict(zip(BrIS_ready, intensity))
    # order detected BrISs by canRT and calculate all sublists
    canRT_ordered_BrIS_ready = [pep for pep, _ in sorted(BrIS_ready_canRT.items(), key=lambda item: item[1])]
    canRT_ordered_BrIS_ready_sublists = calculate_sublists(canRT_ordered_BrIS_ready)
    # find the longest common sublist
    intersections = set(canRT_ordered_BrIS_ready_sublists) & set(refRT_ordered_BrIS_ready_sublists)
    if not len(intersections):
        return (None, None)
    tmp_max_len = max(len(sublist) for sublist in intersections)
    # for all longest sublists, find the candidate with the best correlation
    tmp_best_cor = 0
    tmp_best_candidate = {}
    for sublist in intersections:
        if len(sublist) == tmp_max_len:
            sublist_canRT = np.array([BrIS_ready_canRT[pep] for pep in sublist])
            sublist_refRT = np.array([BrIS_ready_refRT[pep] for pep in sublist])
            sublist_intensity = [BrIS_ready_intensity[pep] for pep in sublist]
            cor = np.corrcoef(sublist_canRT, sublist_refRT)[0, 1]
            if cor > tmp_best_cor:
                tmp_best_cor = cor
                tmp_best_candidate = {
                    'BrIS_sublist': sublist, 
                    'canRT': sublist_canRT.tolist(), 
                    'intensityApex': sublist_intensity, 
                    'rt_cor': cor
                }     
    return (tmp_max_len, tmp_best_candidate) if tmp_best_candidate else (None, None)
    
def verify_BrIS(all_BrISs_envolopes, all_BrISs_refRT, n_procs):
    BrIS_ready = list(all_BrISs_envolopes.keys())
    BrIS_ready_refRT = {pep: all_BrISs_refRT[pep] for pep in BrIS_ready}
    # order detected BrISs by refRT and calculate all sublists
    refRT_ordered_BrIS_ready = [pep for pep, _ in sorted(BrIS_ready_refRT.items(), key=lambda item: item[1])]
    refRT_ordered_BrIS_ready_sublists = calculate_sublists(refRT_ordered_BrIS_ready)
    # The total number of combinatorial possibilities for n detected BrISs, each with six candidates, is 6^n. 
    BrIS_ready_all_canRTs = [[cluster['iso0_rtApex'] for cluster in all_BrISs_envolopes[pep]] for pep in BrIS_ready]
    iter_BrIS_ready_canRTs = itertools.product(*BrIS_ready_all_canRTs)
    BrIS_ready_all_intensities = [[cluster['iso0_intensityApex'] for cluster in all_BrISs_envolopes[pep]] for pep in BrIS_ready]
    iter_BrIS_ready_intensities = itertools.product(*BrIS_ready_all_intensities)
    # create two iterators
    all_candidates = zip(iter_BrIS_ready_canRTs, iter_BrIS_ready_intensities)
    all_candidates_1, all_candidates_2 = itertools.tee(all_candidates, 2)
    # calculate chunk size for imap
    total_tasks = sum(1 for _ in all_candidates_1)
    chunksize = max(1, total_tasks // (4 * n_procs))
    # define partial function, fixing some arguments for multiprocessing
    process_func = partial(
        process_candidates,
        BrIS_ready = BrIS_ready,
        BrIS_ready_refRT = BrIS_ready_refRT,
        refRT_ordered_BrIS_ready_sublists = refRT_ordered_BrIS_ready_sublists
    )
    # use imap to distribute tasks to multiple worker processes
    with Pool(processes = n_procs) as pool:
        results = pool.imap(process_func, all_candidates_2, chunksize = chunksize)
        candidates = [result for result in results if result[0] is not None]
    if not len(candidates):
        return {}
    max_len = max([length for length, _ in candidates])
    candidates = [candidate for length, candidate in candidates if length == max_len]
    verified_BrIS = sorted(candidates, key = lambda x: x['rt_cor'], reverse = True)[0]
    verified_BrIS_to_check = [candidate for candidate in candidates if candidate['BrIS_sublist'] == verified_BrIS['BrIS_sublist']]
    max_total_intensityApex = 0
    best_candidate = {}
    for candidate in verified_BrIS_to_check:
        total_intensityApex = sum(candidate['intensityApex'])
        if total_intensityApex > max_total_intensityApex:
            max_total_intensityApex = total_intensityApex
            best_candidate = candidate
    if best_candidate:
        del best_candidate['rt_cor']
        verified_BrIS = best_candidate
    else:
        verified_BrIS = {}
    return verified_BrIS

# write results
def expand_xic(features, RT_dict, nScans = 5):
    scanApex = features['scanApex']
    scan = features['hill_scan_list']
    intensity = features['hill_intensity_list']
    l = scan[:scan.index(scanApex)]
    r = scan[scan.index(scanApex)+1:]
    if len(l) >= len(r):
        n = len(l) - len(r)
        scan = list(range(scan[0]-nScans, scan[0], 1)) + scan + list(range(scan[-1]+1, scan[-1]+1+n+nScans, 1))
        intensity = [0]*nScans + intensity + [0]*(n+nScans)
    else:
        n = len(r) - len(l)
        scan = list(range(scan[0]-n-nScans, scan[0], 1)) + scan + list(range(scan[-1]+1, scan[-1]+1+nScans, 1))
        intensity = [0]*(n+nScans) + intensity + [0]*nScans
    new_scan, new_rt, new_intensity = [], [], []
    for scan, intv in zip(scan, intensity):
        if scan >= 0 and scan < len(RT_dict):
            new_scan.append(scan)
            new_rt.append(RT_dict[scan])
            new_intensity.append(intv)
    return new_scan, new_rt, new_intensity

def write_results(verified_BrIS, all_BrISs_envolopes, nIsotopes, library_dict, proton_mass, RT_dict, outdir):
    final_results = defaultdict(dict)
    charge_dict = {}
    for pep, rt, intensity in zip(verified_BrIS['BrIS_sublist'], verified_BrIS['canRT'], verified_BrIS['intensityApex']):
        for cluster in all_BrISs_envolopes[pep]:
            if cluster['iso0_rtApex'] == rt and cluster['iso0_intensityApex'] == intensity:
                charge_dict[pep] = cluster['charge']
                for i in range(nIsotopes):
                    iso_tag = f"iso{i}"
                    final_results[iso_tag][pep] = cluster[iso_tag]

    # BrIS all isotopes identification results
    for i in range(nIsotopes):
        iso_tag = f"iso{i}"
        iso_results = final_results[iso_tag]
        basic_anno = {}
        for pep in iso_results:
            sequence = library_dict[pep]['sequence']
            seq_mass = float(library_dict[pep]['isotopes mass'].split(';')[i])
            charge = charge_dict[pep]
            seq_mz = (seq_mass + charge * proton_mass) / charge
            basic_anno[pep] = {'sequence':sequence, 'mass':seq_mass, 'charge':charge, 'mz':seq_mz}
        basic_anno = pd.DataFrame(basic_anno).T
        iso_results = pd.DataFrame(iso_results).T
        iso_results.rename(columns = {'mz': 'hill_mz'}, inplace = True)
        results = basic_anno.join(iso_results)
        results.index.name = 'accession'
        results.to_csv(normalize_path(f"{outdir}/BrIS_isotope{i}_result.csv"))

    # BrIS isotope0 xic
    scans = list(range(len(RT_dict)))
    rts = [rt for rt in RT_dict.values()]
    xics = [['scan'] + scans, ['RT'] + rts, ]
    for pep, features in final_results['iso0'].items():
        scan_list = features['hill_scan_list']
        intensity_list = features['hill_intensity_list']
        scan_intensity = dict(zip(scan_list, intensity_list))
        xic = [scan_intensity.get(idx, 0) for idx in scans]
        xics.append([f"{pep}_intensity"] + xic)        
    with open(normalize_path(f"{outdir}/BrIS_isotope0_xic.csv"), 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',')
        writer.writerows(xics)

    # BrIS all isotopes xis
    all_isotopes_xic = []
    for pep, iso0_features in final_results['iso0'].items():
        iso0_scan, iso0_rt, iso0_intv = expand_xic(iso0_features, RT_dict)
        all_isotopes_xic.extend([[f"{pep}_iso0_rt"] + iso0_rt, [f"{pep}_iso0_intv"] + iso0_intv])
        for i in range(1, nIsotopes):
            iso_tag = f"iso{i}"
            iso_features = final_results[iso_tag][pep]
            iso_scan, iso_rt, iso_intv = expand_xic(iso_features, RT_dict)
            iso_rt_1, iso_intv_1 = [], []
            for scan, rt, intv in zip(iso_scan, iso_rt, iso_intv):
                if scan in iso0_scan:
                    iso_rt_1.append(rt)
                    iso_intv_1.append(intv)
            all_isotopes_xic.extend([[f"{pep}_{iso_tag}_rt"] + iso_rt_1, [f"{pep}_{iso_tag}_intv"] + iso_intv_1])
    with open(normalize_path(f"{outdir}/BrIS_all_isotopes_xic.csv"), 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',')
        writer.writerows(all_isotopes_xic)

    # BrIS all isotopes ratios
    all_isotopes_scan_intv = {}
    for pep in verified_BrIS['BrIS_sublist']:
        local_isotope_scan_intv = {}
        for i in range(nIsotopes):
            iso_tag = f"iso{i}"
            iso_features = final_results[iso_tag][pep]
            local_isotope_scan_intv[iso_tag] = dict(zip(iso_features['hill_scan_list'], iso_features['hill_intensity_list']))
        all_isotopes_scan_intv[pep] = local_isotope_scan_intv

    all_ratios = []
    for pep, iso_features in all_isotopes_scan_intv.items():
        iso_intv_distribution = [float(x) for x in library_dict[pep]['isotopes intensity distribution'].split(';')]
        ratios = {}
        for i in range(1, nIsotopes):
            ratios[f"iso{i}:iso0"] = {'e_ratio':[], 't_ratio':[], 'iso_tag':[]}
        for scan, intv in iso_features['iso0'].items():
            for i in range(1, nIsotopes):
                iso_tag = f"iso{i}"
                if scan in iso_features[iso_tag]:
                    ratios[f"{iso_tag}:iso0"]['e_ratio'].append(iso_features[iso_tag][scan]/intv)
                    ratios[f"{iso_tag}:iso0"]['t_ratio'].append(iso_intv_distribution[i]/iso_intv_distribution[0])
                    ratios[f"{iso_tag}:iso0"]['iso_tag'].append(f"{iso_tag}:iso0")
        ratios_df = [pd.DataFrame(params) for params in ratios.values()]
        ratios_df = pd.concat(ratios_df, ignore_index = True)
        ratios_df['peptide'] = [pep] * len(ratios_df)
        all_ratios.append(ratios_df)
    all_ratios = pd.concat(all_ratios, ignore_index = True)
    all_ratios.to_csv(normalize_path(f"{outdir}/BrIS_isotopes_ratios.csv"), index = False)
    
# ================================ workflow ================================
def run(args, logger):
    # data preparation
    msfile, library, outpath = prepare(args['msfile'], args['library'], args['outpath'])

    # start
    logger.info(f"MS data: {os.path.basename(msfile)}")

    # data format conversion if needed
    if msfile.endswith(".raw"):
        logger.info("Converting .raw to .mzML ...")
        convert_status = convert_raw_to_mzml(msfile)
        if convert_status == 0:
            msfile = f"{os.path.splitext(msfile)[0]}.mzML"
            logger.info(f"Successfully converted to: {msfile}")
        else:
            raise Exception("File format conversion failed!")
        
    # extract ms1
    logger.info("Extracting MS1 ...")
    ms1 = extract_ms1(msfile)

    # process ms1
    logger.info(f"Processing MS1: {len(ms1)} scans ...")
    ms1, max_mz_value = process_ms1(ms1, args['minmz'], args['maxmz'], args['mini'], args['profile'], args['mdiff'])

    # build hills
    mz_step = args['htol'] * max_mz_value * 1e-6
    im_tol = 0 
    logger.info("Building hills ...")
    hills_dict = build_hills(ms1, args['htol'], mz_step, im_tol)

    # split hills
    logger.info(f"Splitting built hills: {len(set(hills_dict['hill_idx_array']))} entries ...")
    hills_dict = split_peaks_multi(hills_dict, ms1, args['hminl'], args['hillValleyFactor'], args['nprocs']) 

    # process hills
    logger.info(f"Processing splitted hills: {len(set(hills_dict['hill_idx_array']))} entries ...")
    hills_dict = process_hills(hills_dict, ms1, mz_step, im_tol)
    RT_dict = {idx:z['retention time'] for idx, z in enumerate(ms1)}
    hills_dict, hills_features = process_hills_extra(hills_dict, RT_dict)

    # detect BrIS envolopes
    library = pd.read_csv(args['library'], index_col = 0)
    library_dict = library.T.to_dict()
    logger.info("Detecting BrIS envolopes ...")
    all_BrISs_envolopes = detect_BrIS_envolopes(library_dict, args['proton_mass'], hills_features, args['ftol'], args['nIsotopes'], args['mincor'], args['nFeatures'])

    # verify detected BrIS envolopes
    all_BrISs_refRT = {BrIS_accession:info['RT'] for BrIS_accession, info in library_dict.items()}
    logger.info("Verifying BrIS ...")
    verified_BrIS = verify_BrIS(all_BrISs_envolopes, all_BrISs_refRT, args['nprocs'])
        
    # write results
    outdir = normalize_path(f"{outpath}/{os.path.basename(msfile).split('.')[0]}")
    os.makedirs(outdir, exist_ok = True)
    logger.info(f"Writing results: {len(verified_BrIS['BrIS_sublist'])} BrIS ...")
    write_results(verified_BrIS, all_BrISs_envolopes, args['nIsotopes'], library_dict, args['proton_mass'], RT_dict, outdir)
    if args['write_hills_features']:
        pd.DataFrame(hills_features).T.to_csv(normalize_path(f"{outdir}/hills_features.csv"), index = False)

    # end
    logger.info("Completed!")

# ================================== main ==================================
if __name__ == "__main__":
    # parameter setting
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(
        description = 'BrIS (Brominated Internal Standard) identification via MS1 spectra.',
        epilog = dedent("""
                        Example usage: 
                            $ cd /PATH/TO/BrIS
                            $ python BrIS_MS1_identification.py -f MSfile.mzML -l ./BrIS_library.csv -o ./outpath
                        Output:
                            - ./outpath/MSfile/BrIS_isotope0_result.csv
                            - ./outpath/MSfile/BrIS_isotope1_result.csv
                            - ./outpath/MSfile/BrIS_isotope2_result.csv
                            - ./outpath/MSfile/BrIS_isotope0_xic.csv
                            - ./outpath/MSfile/BrIS_all_isotopes_xic.csv
                            - ./outpath/MSfile/BrIS_isotopes_ratios.csv
                        """),
        formatter_class = argparse.RawDescriptionHelpFormatter)
    # required parameters
    parser.add_argument('--msfile', 
                        '-f', 
                        help = 'MS file in .d, .raw or .mzML format.', 
                        required = True, 
                        type = str)
    parser.add_argument('--library', 
                        '-l', 
                        help = 'library.csv', 
                        required = True, 
                        type = str)
    parser.add_argument('--outpath', 
                        '-o', 
                        help = 'path to output results', 
                        required = True, 
                        type = str)
    # optional parameters
    parser.add_argument('--minmz', 
                        help = 'minimum m/z in MS1 spectra', 
                        required = False, 
                        default = 350, 
                        type = float)
    parser.add_argument('--maxmz', 
                        help = 'maximum m/z in MS1 spectra', 
                        required = False, 
                        default = 1500, 
                        type = float)
    parser.add_argument('--mini', 
                        help = 'minimum intensity in MS1 spectra', 
                        required = False, 
                        default = 100, 
                        type = float)
    parser.add_argument('--mdiff', 
                        help = 'mass difference (Da) for centroiding peaks', 
                        required = False, 
                        default = 0.05, 
                        type = float)
    parser.add_argument('--htol', 
                        help = 'mass tolerance (ppm) for building hills', 
                        required = False, 
                        default = 8, 
                        type = float)
    parser.add_argument('--ftol', 
                        help = 'mass tolerance (ppm) for detect BrIS envolopes', 
                        required = False, 
                        default = 8, 
                        type = float)
    parser.add_argument('--nIsotopes', 
                        help = 'number of isotopic peaks', 
                        required = False, 
                        default = 3, 
                        type = int)
    parser.add_argument('--nFeatures', 
                        help = 'maximum number of candidate envolopes for each BrIS', 
                        required = False, 
                        default = 3, 
                        type = int)
    parser.add_argument('--mincor', 
                        help = 'minimum correlation coefficient for detecting BrIS envolopes', 
                        required = False, 
                        default = 0.5, 
                        type = float)
    parser.add_argument('--profile', 
                        help = 'whether MS file is in profile mode. if True MS file will be converted into centroided mode.', 
                        action = 'store_true')
    parser.add_argument('--hminl', 
                        help = 'minimum length for hills', 
                        required = False, 
                        default = 2, 
                        type = int)
    parser.add_argument('--hillValleyFactor', 
                        help = 'hill valley factor for splitting hills', 
                        required = False, 
                        default = 1.7, 
                        type = float)
    parser.add_argument('--proton_mass', 
                        help = 'proton mass', 
                        required = False, 
                        default = 1.00727646677, 
                        type = float)
    parser.add_argument('--nprocs', 
                        '-p',
                        help = 'number of process', 
                        required = False, 
                        default = 4, 
                        type = int)
    parser.add_argument('--write_hills_features',
                        help = 'write hills features to csv',
                        action = 'store_true')
    args = vars(parser.parse_args())

    # configure logger
    import logging
    logging.basicConfig(level = logging.INFO, format = "%(asctime)s - %(levelname)s - %(message)s")
    logger = logging.getLogger()

    # run
    run(args, logger)

    
