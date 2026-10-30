#!/usr/bin/env python3
"""Rebuild per-run task events for nm000379 with a Python port of the authors' MATLAB pipeline.

Source code ported (shipped in sourcedata/KragelEtal21_SciAdv/code/, GPL-3.0-or-later, (C) 2021 James Kragel):
  preproc/create_events.m    -> create_events(), get_eye_events(), merge_events()
  preproc/adjust_S1_events.m -> adjust_s1_events()   (participant S1 = sub-01 only)
  utils/get_onsets.m         -> trial_onsets()        (sync-pulse onsets from the DC channel)
  utils/epoch_data.m, bosc/run_bosc_epoch*.m -> alignment rule of eye events to iEEG samples
Inputs (all in the dataset): EyeLink EDF files data/eye/S<N>/S<N>_{s,r}<block>.edf, behavioural arrays
data/behav/S<N>/*, and the DC sync channel (DC03/DC04) of each BIDS run (sub-0N/ieeg/*_ieeg.eeg, the lossless
copy of S<N>_sl_block<k>.m00). EDF files are read with SR Research's EDF Access API (libedfapi.so, the library
the authors' edfmex/Edf2Mat wraps), loaded through ctypes; give its path with --edfapi.

MATLAB semantics kept on purpose: uint32 saturating arithmetic on EyeLink times, round-half-away-from-zero on
uint32 conversion, findpeaks() flat-peak rule (first sample of a plateau) and MinPeakProminence, 1-based indices.
Deviation (documented): the released S<N>_stim_array*.txt files hold COCO file names ("000000067269.jpg"); MATLAB
load() cannot parse them, so they are read as text and converted to the integer COCO id, which is what the authors'
code expects (append_salience.m uses sprintf('%012d', stim_idx)).

Usage: fix-20261030_events.py --root <dataset root with annexed content> --edfapi <libedfapi.so>
                               [--out <dir to write events into, default root>] [--report <json>] [--write]
Without --write only the diagnostics/report are produced.
"""
import argparse
import ctypes
import glob
import json
import math
import os
import re
import sys

import numpy as np

SRC = "sourcedata/KragelEtal21_SciAdv"
SUBJECTS = {1: "S1", 2: "S2", 3: "S3", 4: "S4", 5: "S5", 6: "S6"}
UMAX = 2**32 - 1
N_ENC, N_REC = 24, 48  # create_events.m: trials per phase


# ----------------------------------------------------------------------------------------------------------------
# MATLAB uint32 helpers
def u32(x):
    """MATLAB conversion of a double to uint32: round half away from zero, saturate to [0, 2^32-1]."""
    x = np.asarray(x, dtype=float)
    r = np.sign(x) * np.floor(np.abs(x) + 0.5)
    return np.clip(r, 0, UMAX)


# ----------------------------------------------------------------------------------------------------------------
# EDF reading via libedfapi (same element stream as edfmex: FEVENT codes and RECORDINGS)
CODES = {1: "STARTPARSE", 2: "ENDPARSE", 10: "BREAKPARSE", 3: "STARTBLINK", 4: "ENDBLINK", 5: "STARTSACC",
         6: "ENDSACC", 7: "STARTFIX", 8: "ENDFIX", 9: "FIXUPDATE", 24: "MESSAGEEVENT", 15: "STARTSAMPLES",
         16: "ENDSAMPLES", 17: "STARTEVENTS", 18: "ENDEVENTS"}
SAMPLE_TYPE, RECORDING_INFO, NO_PENDING_ITEMS = 200, 30, 0


class LSTRING(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("len", ctypes.c_int16), ("c", ctypes.c_char)]


_F = ["hstx", "hsty", "gstx", "gsty", "sta", "henx", "heny", "genx", "geny", "ena", "havx", "havy", "gavx",
      "gavy", "ava", "avel", "pvel", "svel", "evel", "supd_x", "eupd_x", "supd_y", "eupd_y"]


class FEVENT(ctypes.Structure):  # edf_data.h, #pragma pack(1) on x86_64
    _pack_ = 1
    _fields_ = ([("time", ctypes.c_uint32), ("type", ctypes.c_int16), ("read", ctypes.c_uint16),
                 ("sttime", ctypes.c_uint32), ("entime", ctypes.c_uint32)]
                + [(n, ctypes.c_float) for n in _F]
                + [("eye", ctypes.c_int16), ("status", ctypes.c_uint16), ("flags", ctypes.c_uint16),
                   ("input", ctypes.c_uint16), ("buttons", ctypes.c_uint16), ("parsedby", ctypes.c_uint16),
                   ("message", ctypes.POINTER(LSTRING))])


class RECORDINGS(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("time", ctypes.c_uint32), ("sample_rate", ctypes.c_float), ("eflags", ctypes.c_uint16),
                ("sflags", ctypes.c_uint16), ("state", ctypes.c_ubyte), ("record_type", ctypes.c_ubyte),
                ("pupil_type", ctypes.c_ubyte), ("recording_mode", ctypes.c_ubyte), ("filter_type", ctypes.c_ubyte),
                ("pos_type", ctypes.c_ubyte), ("eye", ctypes.c_ubyte)]


_LIB = None


def edfapi(path):
    global _LIB
    if _LIB is None:
        lib = ctypes.CDLL(path)
        lib.edf_open_file.restype = ctypes.c_void_p
        lib.edf_open_file.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int)]
        lib.edf_get_next_data.argtypes = [ctypes.c_void_p]
        lib.edf_get_float_data.restype = ctypes.c_void_p
        lib.edf_get_float_data.argtypes = [ctypes.c_void_p]
        lib.edf_get_element_count.argtypes = [ctypes.c_void_p]
        lib.edf_get_element_count.restype = ctypes.c_uint
        lib.edf_close_file.argtypes = [ctypes.c_void_p]
        _LIB = lib
    return _LIB


def read_edf(lib, fname):
    """Return dict(FEVENT=list of dicts in file order, RECORDINGS=list of dicts) like edfmex()."""
    err = ctypes.c_int(0)
    # edfmex defaults: consistency_check=1, load_events=1, load_samples=1
    ef = lib.edf_open_file(fname.encode(), 1, 1, 1, ctypes.byref(err))
    if not ef:
        raise RuntimeError(f"edf_open_file failed for {fname} (err {err.value})")
    n = lib.edf_get_element_count(ef)
    fev, recs = [], []
    for _ in range(n):
        code = lib.edf_get_next_data(ef)
        if code == NO_PENDING_ITEMS:
            break
        if code in CODES:
            p = ctypes.cast(lib.edf_get_float_data(ef), ctypes.POINTER(FEVENT)).contents
            msg = None
            if p.message:
                ln = p.message.contents.len
                addr = ctypes.addressof(p.message.contents) + 2
                msg = ctypes.string_at(addr, max(ln, 0)).split(b"\x00")[0].decode("latin-1")
            d = {"codestring": CODES[code], "sttime": int(p.sttime), "entime": int(p.entime), "message": msg}
            for k in ("gstx", "gsty", "genx", "geny", "gavx", "gavy"):
                d[k] = float(p.__getattribute__(k))
            fev.append(d)
        elif code == RECORDING_INFO:
            r = ctypes.cast(lib.edf_get_float_data(ef), ctypes.POINTER(RECORDINGS)).contents
            recs.append({"time": int(r.time), "state": int(r.state)})
    lib.edf_close_file(ef)
    return {"FEVENT": fev, "RECORDINGS": recs}


# ----------------------------------------------------------------------------------------------------------------
# create_events.m port
def dlmread(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.strip():
                rows.append([float(v) if v.strip() else 0.0 for v in line.rstrip("\n").split("\t")])
    w = max(len(r) for r in rows)
    return np.array([r + [0.0] * (w - len(r)) for r in rows])  # dlmread pads short rows with 0


def load_stim(path):
    rows = []
    with open(path) as f:
        for line in f:
            toks = [t for t in re.split(r"\s+", line.strip()) if t]
            if toks:
                rows.append([int(os.path.splitext(t)[0]) for t in toks])
    return rows  # stim_info(row, col) -> rows[row-1][col-1]


def get_eye_events(edf, block, delay, phase):
    fe = edf["FEVENT"]
    st = np.array([e["sttime"] for e in fe], dtype=float)
    en = np.array([e["entime"] for e in fe], dtype=float)
    code = np.array([e["codestring"] for e in fe])
    msg = np.array([e["message"] if e["message"] is not None else "\x00none" for e in fe])
    start_times = np.array([r["time"] for r in edf["RECORDINGS"] if r["state"] == 1], dtype=float)
    stop_times = np.array([r["time"] for r in edf["RECORDINGS"] if r["state"] == 0], dtype=float)
    out = {"trial": [], "fixation": [], "saccade": [], "blink": []}
    for t in range(len(start_times)):
        lo = u32(start_times[t] + delay[t])  # uint32 + uint32/double -> uint32 (saturating, rounded)
        hi = stop_times[t]
        win = ((en >= lo) & (en < hi)) | ((st >= lo) & (st < hi))
        fix_idx = win & (code == "ENDFIX")
        sac_idx = win & (code == "ENDSACC")
        for i in np.flatnonzero(sac_idx):  # remove ENDSACC produced by blinks
            if i >= 1 and "BLINK" in code[i - 1]:
                sac_idx[i] = False
        blink_idx = win & (code == "ENDBLINK")
        lo_t = u32(start_times[t] - 100)
        win_t = ((en >= lo_t) & (en < hi)) | ((st >= lo_t) & (st < hi))
        trial_idx = win_t & (code == "MESSAGEEVENT") & (msg == phase)
        for name, idx in (("trial", trial_idx), ("fixation", fix_idx), ("saccade", sac_idx), ("blink", blink_idx)):
            ii = np.flatnonzero(idx)
            out[name].append({"phase": phase, "block": block, "trial": t + 1, "type": name,
                              "start_time": st[ii].astype(np.int64).tolist(), "end_time": en[ii].astype(np.int64).tolist(),
                              "start_x": [fe[i]["gstx"] for i in ii], "start_y": [fe[i]["gsty"] for i in ii],
                              "end_x": [fe[i]["genx"] for i in ii], "end_y": [fe[i]["geny"] for i in ii],
                              "mean_x": [fe[i]["gavx"] for i in ii], "mean_y": [fe[i]["gavy"] for i in ii]})
    return out


def create_events(lib, root, subjid):
    data_dir = os.path.join(root, SRC, "data", "behav", subjid) + "/"
    edf_dir = os.path.join(root, SRC, "data", "eye", subjid)
    split_files = not os.path.exists(data_dir + subjid + "_enc_array.txt")
    if not split_files:
        enc_info = dlmread(data_dir + subjid + "_enc_array.txt")
        rec_info = dlmread(data_dir + subjid + "_rec_array.txt")
        stim_info = load_stim(data_dir + subjid + "_stim_array.txt")
    behav = dlmread(data_dir + subjid + "_output_array.txt")
    m = behav[:, 4] != 0
    recog_resp, recog_rt, recog_block = behav[m, 4], behav[m, 1], behav[m, 5]
    warn = []
    if len(recog_resp) % 48 != 0:
        warn.append("Check recognition responses in behavioral output.")
    n_blocks_completed = behav.shape[0] / (N_ENC + N_REC)
    events, diag = [], {}
    for b in range(1, int(math.floor(n_blocks_completed)) + 1):
        if split_files:
            sfx = "_s1.txt" if b < 4 else "_s2.txt"
            enc_info = dlmread(data_dir + subjid + "_enc_array" + sfx)
            rec_info = dlmread(data_dir + subjid + "_rec_array" + sfx)
            stim_info = load_stim(data_dir + subjid + "_stim_array" + sfx)
        edf_s = read_edf(lib, os.path.join(edf_dir, f"{subjid}_s{b}.edf"))
        st = np.array([r["time"] for r in edf_s["RECORDINGS"] if r["state"] == 1], dtype=float)
        sp = np.array([r["time"] for r in edf_s["RECORDINGS"] if r["state"] == 0], dtype=float)
        if len(st) != len(sp):
            raise RuntimeError(f"{subjid} s{b}: start/stop recording count mismatch {len(st)} {len(sp)}")
        delay = u32(u32(sp - st) - 3000)  # uint32 durations - 3000 (saturating)
        enc = get_eye_events(edf_s, b, delay, "encode")
        edf_r = read_edf(lib, os.path.join(edf_dir, f"{subjid}_r{b}.edf"))
        rst = np.array([r["time"] for r in edf_r["RECORDINGS"] if r["state"] == 1], dtype=float)
        rsp = np.array([r["time"] for r in edf_r["RECORDINGS"] if r["state"] == 0], dtype=float)
        rec_duration = recog_rt[recog_block == b]
        if len(rst) != len(rsp) or len(rst) != len(rec_duration):
            raise RuntimeError(f"{subjid} r{b}: recordings {len(rst)}/{len(rsp)} vs responses {len(rec_duration)}")
        rdelay = u32(rsp - rst) - rec_duration  # double(durations) - rt (double, may be negative)
        rec = get_eye_events(edf_r, b, rdelay, "recog")
        diag[b] = {"enc_recordings": int(len(st)), "rec_recordings": int(len(rst)),
                   "enc_msgs": [len(x["start_time"]) for x in enc["trial"]],
                   "rec_msgs": [len(x["start_time"]) for x in rec["trial"]],
                   "enc_rec_start": st.astype(int).tolist(), "enc_rec_stop": sp.astype(int).tolist(),
                   "rec_rec_start": rst.astype(int).tolist(), "rec_rec_stop": rsp.astype(int).tolist(),
                   "msg_texts": sorted({e["message"] for e in edf_s["FEVENT"] + edf_r["FEVENT"]
                                        if e["codestring"] == "MESSAGEEVENT" and e["message"]})[:40]}
        trial = enc["trial"] + rec["trial"]
        fix = enc["fixation"] + rec["fixation"]
        sac = enc["saccade"] + rec["saccade"]
        blink = enc["blink"] + rec["blink"]
        enc_block = np.repeat(np.arange(1, b + 1), N_ENC)
        rec_block = np.repeat(np.arange(1, b + 1), N_REC)
        events = merge_events(events, trial, fix, blink, sac, enc_info, enc_block, rec_info, rec_block,
                              recog_resp, recog_rt, stim_info)
    return events, diag, warn


def merge_events(prev, trial, fix, blink, sac, enc_info, enc_block, rec_info, rec_block, recog_resp, recog_rt,
                 stim_info):
    ev = [e for e in trial + fix + blink + sac if len(e["start_time"]) > 0]
    ev = sorted(ev, key=lambda e: e["start_time"][0])  # MATLAB sort is stable
    prev_stims = {e["stim_idx"] for e in prev}
    for b in sorted({e["block"] for e in ev}):
        enc_out_idx = np.flatnonzero(enc_block == b)  # 0-based rows
        for e in ev:
            if e["block"] != b or e["phase"] != "encode":
                continue
            idx = enc_out_idx[e["trial"] - 1]
            row = enc_info[idx]
            e["stim_idx"] = stim_info[int(row[1]) - 1][int(row[0]) - 1]
            e["condition"] = int(row[1])
            old = int(row[2])
            if prev:
                old = int(bool(old) or (e["stim_idx"] in prev_stims))
            e["old"] = old
            e["salience"] = int(row[3])
            e["resp"] = float("nan")
            e["rt"] = float("nan")
        rec_out_idx = np.flatnonzero(rec_block == b)
        for e in ev:
            if e["block"] != b or e["phase"] != "recog":
                continue
            idx = rec_out_idx[e["trial"] - 1]
            row = rec_info[idx]
            e["stim_idx"] = stim_info[int(row[1]) - 1][int(row[0]) - 1]
            e["condition"] = int(row[1])
            old = int(row[2])
            if prev:
                old = int(bool(old) or (e["stim_idx"] in prev_stims))
            e["old"] = old
            e["salience"] = int(row[3])
            e["resp"] = float(recog_resp[idx])
            e["rt"] = float(recog_rt[idx])
    return ev + prev


# ----------------------------------------------------------------------------------------------------------------
# adjust_S1_events.m port (MATLAB integer-class behaviour reproduced)
def adjust_s1_events(events):
    """Returns (events with adjusted recog start times, adjust matrix [n_block x 48] as MATLAB builds it)."""
    blocks = sorted({e["block"] for e in events})
    rows = []
    out = []
    for b in blocks:
        bev = [dict(e) for e in events if e["block"] == b]
        ids = sorted({e["stim_idx"] for e in bev})
        adjust = {}
        for i, sid in enumerate(ids, start=1):
            tr = [e for e in bev if e["phase"] == "recog" and e["stim_idx"] == sid and e["type"] == "trial"]
            if not tr:
                continue
            tr = tr[0]
            if len(tr["start_time"]) != 1:
                raise RuntimeError(f"S1 block {b} stim {sid}: {len(tr['start_time'])} recog messages")
            # uint32(end) - uint32(start) saturates at 0; minus double (rt + 800) -> uint32, saturating at 0
            d = max(tr["end_time"][0] - tr["start_time"][0], 0)
            a = int(u32(d - (tr["rt"] + 800)))
            adjust[i] = a
            tr["start_time"] = [int(u32(tr["start_time"][0] + a))]
            for typ in ("fixation", "saccade", "blink"):
                for e in bev:
                    if e["phase"] == "recog" and e["stim_idx"] == sid and e["type"] == typ:
                        e["start_time"] = [int(u32(x - a)) for x in e["start_time"]]
        n = max(adjust) if adjust else 0
        rows.append([adjust.get(i, 0) for i in range(1, n + 1)])  # uint32 array grown by index, gaps = 0
        out += bev
    width = max(48, max(len(r) for r in rows))
    mat = np.zeros((len(rows), width))
    for k, r in enumerate(rows):
        mat[k, :len(r)] = r
    return out, mat


# ----------------------------------------------------------------------------------------------------------------
# get_onsets.m port: findpeaks(sync, 'MinPeakProminence', 2); onsets = locs - 1
def matlab_findpeaks_prominence(y, min_prom):
    y = np.asarray(y, dtype=float)
    n = len(y)
    # keep first of each run of equal values (flat peaks -> lowest index), then sign of derivative
    keep = np.r_[True, y[1:] != y[:-1]]
    iv = np.flatnonzero(keep)
    yv = y[iv]
    # padded with NaN at both ends in MATLAB -> endpoints are never peaks
    s = np.sign(np.diff(yv))
    imax = np.flatnonzero(np.diff(s) < 0) + 1
    pk = iv[imax]
    if len(pk) == 0:
        return pk
    # prominence: extend left/right until a strictly higher sample or the signal end; base = max of the two minima
    prom = np.empty(len(pk))
    for k, p in enumerate(pk):
        h = y[p]
        j = p - 1
        lmin = h
        while j >= 0 and y[j] <= h:
            lmin = min(lmin, y[j]); j -= 1
        j = p + 1
        rmin = h
        while j < n and y[j] <= h:
            rmin = min(rmin, y[j]); j += 1
        prom[k] = h - max(lmin, rmin)
    return pk[prom >= min_prom]


def trial_onsets(sync_counts):
    """sync_counts: DC channel in BrainVision INT_32 counts (0.01 uV); returns 1-based MATLAB sample indices."""
    # MATLAB sees the ASCII microvolt values (= counts / 100); prominence 2 -> 200 counts, exact in integers
    pk0 = matlab_findpeaks_prominence(sync_counts, 200)
    locs = pk0 + 1  # MATLAB 1-based
    onsets = locs - 1
    enc_idx = [i for i in range(3, 50, 2) if i <= len(onsets)]
    rec_idx = [i for i in range(54, 102) if i <= len(onsets)]
    return onsets, np.array([onsets[i - 1] for i in enc_idx]), np.array([onsets[i - 1] for i in rec_idx])


def read_dc(vhdr):
    txt = open(vhdr, encoding="utf-8").read()
    nch = int(re.search(r"NumberOfChannels=(\d+)", txt).group(1))
    assert "BinaryFormat=INT_32" in txt and "DataOrientation=MULTIPLEXED" in txt
    chans = re.findall(r"^Ch(\d+)=([^,]*),", txt, flags=re.M)
    dc = [int(i) for i, name in chans if name.startswith("DC") or name.startswith("*DC")]
    if len(dc) != 1:  # ecog_reref.m: error('Could not find correct sync channel')
        raise RuntimeError(f"{vhdr}: {len(dc)} DC channels")
    fs = 1e6 / float(re.search(r"SamplingInterval=([\d.]+)", txt).group(1))
    eeg = os.path.join(os.path.dirname(vhdr), re.search(r"DataFile=(.*)", txt).group(1).strip())
    data = np.memmap(eeg, dtype="<i4", mode="r")
    assert data.size % nch == 0
    return np.array(data.reshape(-1, nch)[:, dc[0] - 1]), fs, data.size // nch


def interval_check(sync_onsets, msg_by_trial, fs):
    """Residuals (ms) between sync-pulse onset intervals and EyeLink trial-message intervals, paired by trial
    number (onset k <-> trial k, the authors' rule), both relative to the first paired trial."""
    pairs = [(t, msg_by_trial[t]) for t in sorted(msg_by_trial) if t <= len(sync_onsets)]
    if len(pairs) < 2:
        return None
    s = np.array([sync_onsets[t - 1] for t, _ in pairs], float) / fs * 1000.0
    m = np.array([v for _, v in pairs], float)
    r = (s - s[0]) - (m - m[0])
    return {"n": len(pairs), "max_abs_ms": float(np.max(np.abs(r))), "resid_ms": [round(float(x), 2) for x in r]}


def msgs(events, blk, phase):
    return {e["trial"]: e["start_time"][0] for e in events
            if e["block"] == blk and e["type"] == "trial" and e["phase"] == phase and len(e["start_time"]) == 1}


# ----------------------------------------------------------------------------------------------------------------
# Assembly into BIDS events.tsv (alignment rule of utils/epoch_data.m and bosc/run_bosc_epoch*.m):
#   trial onset (samples, 1-based)  encode: enc_onsets(block, trial); recog: rec_onsets(block, trial) + rec_adjust/1000*fs
#   eye event onset (samples)       ref_idx + (eye.start_time - trial.start_time)/1000*fs, ref_idx = un-adjusted onset
# BIDS onset (s) = (1-based sample - 1) / fs; sample column is 0-based.
MAX_RESID_MS = 50.0  # acceptance: sync-pulse vs EyeLink-message interval residuals (pulse/trial mapping check)
COLS = ["onset", "duration", "trial_type", "sample", "phase", "block", "trial", "stim_id", "condition", "old",
        "cue_salience", "response", "response_time", "eyelink_time", "start_x", "start_y", "end_x", "end_y",
        "mean_x", "mean_y"]


def fmt(v, nd=None):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "n/a"
    if nd is not None:
        return f"{v:.{nd}f}"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def build_rows(events, blk, enc_on, rec_on, rec_adjust_row, fs, phases):
    rows = []
    trials = {(e["phase"], e["trial"]): e for e in events if e["block"] == blk and e["type"] == "trial"}
    for e in events:
        if e["block"] != blk or e["phase"] not in phases:
            continue
        key = (e["phase"], e["trial"])
        on = enc_on if e["phase"] == "encode" else rec_on
        if e["trial"] > len(on) or key not in trials:  # run stopped early / no trial message
            continue
        ref_idx = float(on[e["trial"] - 1])
        tstart = trials[key]["start_time"][0]
        common = {"phase": e["phase"], "block": blk, "trial": e["trial"], "stim_id": e["stim_idx"],
                  "condition": e["condition"], "old": e["old"], "cue_salience": e["salience"],
                  "response": e["resp"], "response_time": e["rt"] / 1000.0}  # ms -> s (BIDS unit)
        if e["type"] == "trial":
            adj = 0.0
            if e["phase"] == "recog" and rec_adjust_row is not None:
                adj = float(rec_adjust_row[e["trial"] - 1]) / 1000 * fs
            s1 = ref_idx + adj
            rows.append(dict(common, onset=(s1 - 1) / fs, duration=float("nan"), trial_type=f"{e['phase']}/trial",
                             sample=s1 - 1, eyelink_time=tstart))
            continue
        for k in range(len(e["start_time"])):
            s1 = ref_idx + (e["start_time"][k] - tstart) / 1000 * fs
            rows.append(dict(common, onset=(s1 - 1) / fs, duration=(e["end_time"][k] - e["start_time"][k]) / 1000,
                             trial_type=f"{e['phase']}/{e['type']}", sample=s1 - 1,
                             eyelink_time=e["start_time"][k], start_x=e["start_x"][k], start_y=e["start_y"][k],
                             end_x=e["end_x"][k], end_y=e["end_y"][k], mean_x=e["mean_x"][k],
                             mean_y=e["mean_y"][k]))
            # EyeLink stores no average gaze for saccades and no gaze positions for blinks (0 in every EDF
            # record of this release): written as n/a instead of 0
            if e["type"] == "saccade":
                rows[-1]["mean_x"] = rows[-1]["mean_y"] = None
            elif e["type"] == "blink":
                for c in ("start_x", "start_y", "end_x", "end_y", "mean_x", "mean_y"):
                    rows[-1][c] = None
    order = {"trial": 0, "fixation": 1, "saccade": 2, "blink": 3}
    rows.sort(key=lambda r: (r["onset"], order[r["trial_type"].split("/")[1]]))
    return rows


def write_tsv(path, rows):
    with open(path, "w", newline="\n") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            out = []
            for c in COLS:
                v = r.get(c)
                if c in ("onset", "duration"):
                    out.append(fmt(v, 4))
                elif c == "response_time":
                    out.append(fmt(v, 3))
                elif c == "sample":
                    out.append(str(int(round(v))))
                elif c in ("start_x", "start_y", "end_x", "end_y", "mean_x", "mean_y"):
                    out.append("n/a" if v is None else f"{v:.1f}")
                else:
                    out.append(fmt(v))
            f.write("\t".join(out) + "\n")


def run(a):
    lib = edfapi(a.edfapi)
    report = {}
    for sn in [int(x) for x in a.subjects.split(",")]:
        subj = SUBJECTS[sn]
        events, diag, warn = create_events(lib, a.root, subj)
        rep = {"warnings": warn, "edf": {str(k): v for k, v in diag.items()}, "runs": {}}
        adj = None
        if subj == "S1":  # get_onsets.m: rec_adjust from adjust_S1_events for S1, zeros otherwise
            _, adj = adjust_s1_events(events)
            rep["s1_adjust"] = adj.tolist()
        for vhdr in sorted(glob.glob(os.path.join(a.root, f"sub-{sn:02d}", "ieeg", "*_ieeg.vhdr"))):
            blk = int(re.search(r"run-(\d+)", vhdr).group(1))
            cache = os.path.join(a.dc_cache, os.path.basename(vhdr) + ".dc.npz") if a.dc_cache else None
            if cache and os.path.exists(cache):
                z = np.load(cache); dc, fs, nsamp = z["dc"], float(z["fs"]), int(z["nsamp"])
            else:
                dc, fs, nsamp = read_dc(vhdr)
                if cache:
                    os.makedirs(a.dc_cache, exist_ok=True); np.savez(cache, dc=dc, fs=fs, nsamp=nsamp)
            onsets, enc_on, rec_on = trial_onsets(dc)
            enc_msg, rec_msg = msgs(events, blk, "encode"), msgs(events, blk, "recog")
            cross = {}
            for ob in sorted(diag):  # diagnostic only: how well this run's pulses fit each EyeLink block
                ce = interval_check(enc_on, msgs(events, ob, "encode"), fs)
                cr = interval_check(rec_on, msgs(events, ob, "recog"), fs)
                cross[str(ob)] = [ce and round(ce["max_abs_ms"], 1), cr and round(cr["max_abs_ms"], 1)]
            scan = {}
            for ph, step, n, mm in (("enc", 2, N_ENC, enc_msg), ("rec", 1, N_REC, rec_msg)):
                best = None  # diagnostic only: best first-pulse index if the fixed index rule were shifted
                for k in range(1, len(onsets) + 1):
                    idx = [k + step * i for i in range(n) if k + step * i <= len(onsets)]
                    c = interval_check(np.array([onsets[i - 1] for i in idx]), mm, fs)
                    if c and c["n"] >= 10 and (best is None or c["max_abs_ms"] < best[1]):
                        best = (k, round(c["max_abs_ms"], 1), c["n"])
                scan[ph] = best
            rep["runs"][str(blk)] = {"shift_scan_first_pulse_maxres_n": scan, "fs": fs, "n_samples": int(nsamp), "n_pulses": int(len(onsets)),
                                     "pulse_onsets": onsets.tolist(),
                                     "n_enc_onsets": int(len(enc_on)), "n_rec_onsets": int(len(rec_on)),
                                     "enc_check": interval_check(enc_on, enc_msg, fs),
                                     "rec_check": interval_check(rec_on, rec_msg, fs),
                                     "cross_block_max_abs_ms": cross}
            rr = rep["runs"][str(blk)]
            ok, why = True, []
            if subj == "S6" and blk == 4:  # get_onsets.m: 'missing encoding data for half of block 4, nan onsets'
                ok = False; why.append("authors set S6 block 4 onsets to NaN (get_onsets.m)")
            if blk not in diag:
                ok = False; why.append("no behavioural/EyeLink block")
            else:
                if diag[blk]["enc_recordings"] != N_ENC or diag[blk]["rec_recordings"] != N_REC:
                    ok = False; why.append("EyeLink recording count != 24/48")
                if any(x > 1 for x in diag[blk]["enc_msgs"] + diag[blk]["rec_msgs"]):
                    ok = False; why.append("trial with more than one EyeLink trial message")
                miss = [f"encode {t+1}" for t, x in enumerate(diag[blk]["enc_msgs"]) if x == 0] + \
                       [f"recog {t+1}" for t, x in enumerate(diag[blk]["rec_msgs"]) if x == 0]
                if miss:
                    why.append("no EyeLink trial message (trial and its eye events not aligned, as in the authors' code): " + ", ".join(miss))
            if len(onsets) not in (101, 102):
                ok = False; why.append(f"{len(onsets)} sync pulses detected (authors' index rule assumes 102: 3:2:49 encode, 54:101 recog)")
            phase_ok = {}
            for ph, n_on, n_tr in (("enc", len(enc_on), N_ENC), ("rec", len(rec_on), N_REC)):
                chk = rr[f"{ph}_check"]
                phase_ok[ph] = True
                if n_on != n_tr:
                    phase_ok[ph] = False; why.append(f"{ph}: {n_on} sync onsets vs {n_tr} trials")
                if chk is None or chk["max_abs_ms"] > MAX_RESID_MS:
                    phase_ok[ph] = False; why.append(f"{ph}: interval residual check failed")
            # The authors' rule treats each sync onset as simultaneous with the EyeLink trial message (eye events are
            # placed at onset + (eye time - message time)). Check that the pulse-minus-message offset is the same at
            # recognition as at encoding (same clock mapping within the run).
            ofs = {}
            for ph, on, mm in (("enc", enc_on, enc_msg), ("rec", rec_on, rec_msg)):
                d = [on[t - 1] / fs * 1000.0 - v for t, v in mm.items() if t <= len(on)]
                ofs[ph] = float(np.median(d)) if d else None
            rr["rec_minus_enc_offset_ms"] = (None if None in ofs.values() else round(ofs["rec"] - ofs["enc"], 1))
            if rr["rec_minus_enc_offset_ms"] is None or abs(rr["rec_minus_enc_offset_ms"]) > MAX_RESID_MS:
                phase_ok["rec"] = False
                why.append(f"rec: sync onsets are offset from the EyeLink 'recog' messages by "
                           f"{rr['rec_minus_enc_offset_ms']} ms relative to encoding (pulse timing differs at retrieval)")
            if not phase_ok["enc"]:
                ok = False
            rr["accepted"], rr["phases_written"], rr["reasons"] = ok, [p for p in ("encode", "recog") if ok and phase_ok[p[:3]]], why
            if ok and a.write:
                arow = adj[blk - 1] if adj is not None else None
                rows = build_rows(events, blk, enc_on, rec_on, arow, fs, set(rr["phases_written"]))
                outdir = os.path.join(a.out or a.root, f"sub-{sn:02d}", "ieeg")
                os.makedirs(outdir, exist_ok=True)
                p = os.path.join(outdir, os.path.basename(vhdr).replace("_ieeg.vhdr", "_events.tsv"))
                write_tsv(p, rows)
                rr["written"] = os.path.relpath(p, a.out or a.root)
                rr["n_rows"] = {t: sum(1 for r in rows if r["trial_type"] == t) for t in sorted({r["trial_type"] for r in rows})}
                rr["min_onset_s"] = min(r["onset"] for r in rows)
                rr["max_onset_s"] = max(r["onset"] for r in rows)
                rr["duration_s"] = nsamp / fs
            print(subj, blk, len(onsets), rep["runs"][str(blk)]["enc_check"] and rep["runs"][str(blk)]["enc_check"]["max_abs_ms"],
                  rep["runs"][str(blk)]["rec_check"] and rep["runs"][str(blk)]["rec_check"]["max_abs_ms"], flush=True)
        report[subj] = rep
        if a.dump:
            os.makedirs(a.dump, exist_ok=True)
            with open(os.path.join(a.dump, f"{subj}_events_eyelinktime.json"), "w") as f:
                json.dump(events, f)
    with open(a.report, "w") as f:
        json.dump(report, f, indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--edfapi", required=True)
    ap.add_argument("--out")
    ap.add_argument("--report", required=True)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--dc-cache", help="optional directory caching the DC channel of each run (speed only)")
    ap.add_argument("--dump", help="directory for the per-subject event structures in EyeLink time (JSON)")
    ap.add_argument("--subjects", default="1,2,3,4,5,6")
    a = ap.parse_args()
    run(a)
