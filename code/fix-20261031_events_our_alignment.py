#!/usr/bin/env python3
"""Recover task events for the nm000379 runs that the authors' fixed pulse indices cannot align (OUR alignment).

THIS IS OUR ALIGNMENT, NOT THE AUTHORS'. The main per-run *_events.tsv files (code/fix-20261030_events.py) follow the
authors' rule (sync pulse k = trial k, fixed indices 3:2:49 / 54:101 from utils/get_onsets.m). For runs where that
rule does not apply, this script estimates the mapping between the EyeLink clock and the iEEG run clock itself and
writes the result as a derivative (derivatives/events-our-alignment/), never into the main tables.

Method
  * Inputs and event content are those of the main pipeline: EyeLink EDF + behavioural arrays parsed by the port of
    the authors' create_events.m (imported from code/fix-20261030_events.py), sync pulses found as in get_onsets.m
    (findpeaks, prominence 2 uV, one sample before the peak) on the run's DC channel.
  * Pairing by search (no fixed indices): for one EyeLink block and phase ('encode' or 'recog' trial messages, EyeLink
    ms) and the run's pulses (s), every (message i, pulse j) anchor defines an offset; each message is paired with the
    nearest pulse within TOL (50 ms); the anchor with most pairs wins (tie: lowest RMS). Then a least-squares linear
    clock fit t_run = a + b * t_eyelink is refitted and the messages re-paired twice.
  * Every event is placed through the fitted map: onset = a + b * eyelink_time (trial rows: the EyeLink trial message;
    eye events: their EyeLink start time). Durations are EyeLink durations.
  * sub-01 recognition phase: its pulses lag the 'recog' messages by 0.8-1.2 s (other participants: simultaneous), so
    they are not used; the map fitted on the same run's encoding phase is applied to the recognition messages.

Acceptance criteria (stated before the left-out runs were examined; see the derivative README)
  C1 every message of the phase whose mapped time lies inside the recording is paired, and at least MIN_PAIRS pairs
  C2 max |residual| <= TOL after the fit
  C3 unique: the runner-up (any other EyeLink block of the participant, same phase; or the same block with a pulse
     assignment whose offset differs by > 0.1 s) pairs at most RUNNER_UP_FRAC of the best pair count
  C4 |b - 1| <= MAX_PPM * 1e-6
  C5 if both phases of a run are paired, their offsets agree within TOL and one joint fit is used
  sub-01 recognition: C1-C4 on encoding, plus S1a (recognition pulses minus mapped messages: spread <= TOL) and S1b
  (the encoding-only extrapolation reproduces the recognition pulses of the other participants' accepted runs within
  TOL).
Calibration: the same search is run on the runs the main tables cover; it must recover the authors' indices and agree
with the main onsets (reported in the receipt).

Usage: fix-20261031_events_our_alignment.py --root <dataset root with sourcedata + sub-*/ieeg/*.vhdr (+ .eeg unless
       --dc-cache holds the DC channel)> --edfapi <libedfapi.so> --out <derivative dir> --receipt <json>
       [--main-root <tree with the main events.tsv, for the agreement check>] [--dc-cache <dir>]
"""
import argparse
import glob
import importlib.util
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("port", os.path.join(HERE, "fix-20261030_events.py"))
port = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(port)

TOL = 0.050          # s, pairing tolerance and max residual (C1, C2, C5, S1a, S1b)
MIN_PAIRS = 10       # C1
RUNNER_UP_FRAC = 0.5  # C3
MAX_PPM = 1000.0     # C4
DIFF_ASSIGN = 0.100  # s, offsets differing by more than this are a different pulse assignment (C3)
PRECUE_GAP = (0.70, 0.85)  # s, encoding pre-cue pulse -> image pulse (get_onsets.m: ~750 ms; measured 766-769 ms)
SUBJECTS = {1: "S1", 2: "S2", 3: "S3", 4: "S4", 5: "S5", 6: "S6"}
# runs the main tables cover (lane B, code/fix-20261030_events.py): (participant, run) -> phases
MAIN_RUNS = {**{(1, r): ("encode",) for r in range(1, 9)}, **{(2, r): ("encode", "recog") for r in range(1, 9)},
             **{(4, r): ("encode", "recog") for r in (1, 2, 3)}, **{(5, r): ("encode", "recog") for r in (2, 3, 4)},
             (6, 3): ("encode", "recog"), (6, 5): ("encode", "recog")}
# runs/phases left out of the main tables, examined here
TARGETS = {(1, r): ("recog",) for r in range(1, 9)}
TARGETS.update({(3, 1): ("encode", "recog"), (3, 2): ("encode", "recog"), (3, 4): ("encode", "recog"),
                (5, 1): ("encode", "recog"), (6, 1): ("encode", "recog"), (6, 2): ("encode", "recog"),
                (6, 4): ("encode", "recog")})
COLS = port.COLS + ["alignment"]


def pair(msg_s, pulses, a, b, tol=TOL):
    """Pair each message (EyeLink s) with the nearest pulse (run s) within tol under t = a + b*m."""
    pred = a + b * msg_s
    j = np.clip(np.searchsorted(pulses, pred), 1, len(pulses) - 1)
    jj = np.where(np.abs(pulses[j - 1] - pred) <= np.abs(pulses[j] - pred), j - 1, j)
    d = pulses[jj] - pred
    ok = np.abs(d) <= tol
    return jj, d, ok


def search(msg_s, pulses, b=1.0):
    """All anchors: returns array of (offset, n_pairs, rms) for every (message, pulse) anchor."""
    offs = (pulses[None, :] - b * msg_s[:, None]).ravel()
    offs = np.unique(np.round(offs, 4))
    pred = offs[:, None] + b * msg_s[None, :]
    j = np.clip(np.searchsorted(pulses, pred), 1, len(pulses) - 1)
    d = np.minimum(np.abs(pulses[j - 1] - pred), np.abs(pulses[j] - pred))
    ok = d <= TOL
    n = ok.sum(1)
    rms = np.sqrt(np.where(ok, d ** 2, 0).sum(1) / np.maximum(n, 1))
    return offs, n, rms


def fit(msg_s, pulses, a0, b0=1.0, iters=2):
    a, b = a0, b0
    for _ in range(iters + 1):
        jj, d, ok = pair(msg_s, pulses, a, b)
        if ok.sum() < 2:
            break
        x, y = msg_s[ok], pulses[jj[ok]]
        if ok.sum() >= 3 and np.ptp(x) > 10:
            b, a = np.polyfit(x, y, 1)
        else:
            a = float(np.mean(y - b * x))
    jj, d, ok = pair(msg_s, pulses, a, b)
    return a, b, jj, d, ok


def best_candidate(msg_s, pulses):
    offs, n, rms = search(msg_s, pulses)
    if len(offs) == 0:
        return None
    order = np.lexsort((rms, -n))
    return offs, n, rms, order


def precue_mask(pulses):
    """Pulses followed by another pulse PRECUE_GAP later (encoding pre-cue -> image pulse, get_onsets.m)."""
    lo = np.searchsorted(pulses, pulses + PRECUE_GAP[0], side="left")
    hi = np.searchsorted(pulses, pulses + PRECUE_GAP[1], side="right")
    return hi > lo


def phase_fit(msg_by_trial, pulses_all, dur, phase):
    """Fit one block/phase. Returns dict with fit, pairing and the anchor table for runner-up evaluation.
    Encoding: only pre-cue pulses (followed by the image pulse 0.70-0.85 s later) are pairing targets, because each
    study trial has two pulses and a one-pulse shift would otherwise pair every message as well."""
    sub = np.flatnonzero(precue_mask(pulses_all)) if phase == "encode" else np.arange(len(pulses_all))
    if len(sub) < 2:
        return None
    pulses = pulses_all[sub]
    trials = sorted(msg_by_trial)
    m = np.array([msg_by_trial[t] for t in trials], float) / 1000.0
    bc = best_candidate(m, pulses)
    if bc is None:
        return None
    offs, n, rms, order = bc
    k = order[0]
    a, b, jj, d, ok = fit(m, pulses, offs[k])
    jj = sub[jj]  # indices into all pulses
    pred = a + b * m
    inside = (pred >= 0) & (pred < dur)
    res = d[ok]
    return {"trials": trials, "msg_s": m, "a": float(a), "b": float(b), "pulse_idx": jj, "resid": d, "paired": ok,
            "inside": inside, "n_pairs": int(ok.sum()), "n_msgs": len(trials), "n_inside": int(inside.sum()),
            "max_abs_ms": float(np.max(np.abs(res)) * 1000) if len(res) else None,
            "rms_ms": float(np.sqrt(np.mean(res ** 2)) * 1000) if len(res) else None,
            "ppm": float((b - 1) * 1e6), "anchor_n": n, "anchor_offs": offs, "best_anchor_off": float(offs[k])}


def runner_up(cands, best_key, best_off):
    """Largest pair count among other blocks (same phase) or the same block at a different offset."""
    best = 0
    where = None
    for key, c in cands.items():
        if c is None:
            continue
        if key == best_key:
            sel = np.abs(c["anchor_offs"] - best_off) > DIFF_ASSIGN
        else:
            sel = np.ones(len(c["anchor_offs"]), bool)
        if sel.any():
            i = int(np.argmax(np.where(sel, c["anchor_n"], -1)))
            if c["anchor_n"][i] > best:
                best, where = int(c["anchor_n"][i]), {"block": key, "offset_s": round(float(c["anchor_offs"][i]), 3)}
    return best, where


def load_pulses(vhdr, dc_cache):
    cache = os.path.join(dc_cache, os.path.basename(vhdr) + ".dc.npz") if dc_cache else None
    if cache and os.path.exists(cache):
        z = np.load(cache)
        dc, fs, nsamp = z["dc"], float(z["fs"]), int(z["nsamp"])
    else:
        dc, fs, nsamp = port.read_dc(vhdr)
    onsets, enc_on, rec_on = port.trial_onsets(dc)
    return (onsets.astype(float) - 1) / fs, fs, nsamp, dc  # onset (s) = (1-based sample - 1)/fs, as in the main tables


def build_rows(events, blk, phases, maps, fs, dur, label, pulse_rows=None):
    """Rows for block blk, phases -> map (a, b). Trial rows at the mapped message; eye events at mapped start.
    pulse_rows (sub-01 recognition): {trial: paired sync-pulse time (s)} -> extra 'recog/sync_pulse' rows."""
    rows, dropped = [], 0
    trials = {(e["phase"], e["trial"]): e for e in events if e["block"] == blk and e["type"] == "trial"}
    for e in events:
        if e["block"] != blk or e["phase"] not in phases or (e["phase"], e["trial"]) not in trials:
            continue
        a, b = maps[e["phase"]]
        tmsg = trials[(e["phase"], e["trial"])]["start_time"][0]
        if not (0 <= a + b * tmsg / 1000.0 < dur):  # trial outside the recording (partial runs)
            dropped += 1 if e["type"] == "trial" else 0
            continue
        common = {"phase": e["phase"], "block": blk, "trial": e["trial"], "stim_id": e["stim_idx"],
                  "condition": e["condition"], "old": e["old"], "cue_salience": e["salience"],
                  "response": e["resp"], "response_time": e["rt"] / 1000.0, "alignment": label[e["phase"]]}
        if e["type"] == "trial":
            t = a + b * tmsg / 1000.0
            rows.append(dict(common, onset=t, duration=float("nan"), trial_type=f"{e['phase']}/trial",
                             sample=t * fs, eyelink_time=tmsg))
            if pulse_rows and e["phase"] == "recog" and e["trial"] in pulse_rows:
                rows.append(dict(common, onset=pulse_rows[e["trial"]], duration=float("nan"),
                                 trial_type="recog/sync_pulse", sample=0, eyelink_time=None,
                                 alignment="sync_pulse_paired"))
            continue
        for k in range(len(e["start_time"])):
            t = a + b * e["start_time"][k] / 1000.0
            if not (0 <= t < dur):
                continue
            r = dict(common, onset=t, duration=(e["end_time"][k] - e["start_time"][k]) / 1000,
                     trial_type=f"{e['phase']}/{e['type']}", sample=t * fs, eyelink_time=e["start_time"][k],
                     start_x=e["start_x"][k], start_y=e["start_y"][k], end_x=e["end_x"][k], end_y=e["end_y"][k],
                     mean_x=e["mean_x"][k], mean_y=e["mean_y"][k])
            if e["type"] == "saccade":
                r["mean_x"] = r["mean_y"] = None
            elif e["type"] == "blink":
                for c in ("start_x", "start_y", "end_x", "end_y", "mean_x", "mean_y"):
                    r[c] = None
            rows.append(r)
    order = {"trial": 0, "sync_pulse": 1, "fixation": 2, "saccade": 3, "blink": 4}
    for r in rows:  # write onset and sample consistently: sample = round(onset*fs), onset = sample/fs
        r["sample"] = int(round(r["onset"] * fs))
        r["onset"] = r["sample"] / fs
    rows.sort(key=lambda r: (r["onset"], order[r["trial_type"].split("/")[1]]))
    return rows, dropped


def write_tsv(path, rows):
    with open(path, "w", newline="\n") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            out = []
            for c in COLS:
                v = r.get(c)
                if c in ("onset", "duration"):
                    out.append(port.fmt(v, 4))
                elif c == "response_time":
                    out.append(port.fmt(v, 3))
                elif c == "sample":
                    out.append(str(int(v)))
                elif c in ("start_x", "start_y", "end_x", "end_y", "mean_x", "mean_y"):
                    out.append("n/a" if v is None else f"{v:.1f}")
                else:
                    out.append(port.fmt(v))
            f.write("\t".join(out) + "\n")


def r1(x, n=1):
    return None if x is None else round(float(x), n)


def main(a):
    lib = port.edfapi(a.edfapi)
    receipt = {"criteria": {"tol_ms": TOL * 1000, "min_pairs": MIN_PAIRS, "runner_up_frac": RUNNER_UP_FRAC,
                            "max_ppm": MAX_PPM, "diff_assign_ms": DIFF_ASSIGN * 1000}, "runs": {}, "s1b": []}
    s1_checks = []
    for sn, subj in SUBJECTS.items():
        events, diag, _ = port.create_events(lib, a.root, subj)
        blocks = sorted(diag)
        msgs = {(b, ph): port.msgs(events, b, ph) for b in blocks for ph in ("encode", "recog")}
        for vhdr in sorted(glob.glob(os.path.join(a.root, f"sub-{sn:02d}", "ieeg", "*_ieeg.vhdr"))):
            run = int(re.search(r"run-(\d+)", vhdr).group(1))
            if (sn, run) not in MAIN_RUNS and (sn, run) not in TARGETS:
                continue
            pulses, fs, nsamp, dc = load_pulses(vhdr, a.dc_cache)
            dur = nsamp / fs
            rep = {"participant": f"sub-{sn:02d}", "run": run, "n_pulses": len(pulses), "duration_s": dur,
                   "role": "calibration (main table)" if (sn, run) in MAIN_RUNS else "target", "phases": {}}
            fits = {}
            for ph in ("encode", "recog"):
                cands = {b: (phase_fit(msgs[(b, ph)], pulses, dur, ph) if len(msgs[(b, ph)]) >= 2 else None)
                         for b in blocks}
                ok_c = {b: c for b, c in cands.items() if c is not None}
                if not ok_c:
                    continue
                bb = max(ok_c, key=lambda b: (ok_c[b]["n_pairs"], -(ok_c[b]["rms_ms"] or 1e9)))
                c = ok_c[bb]
                ru_n, ru_where = runner_up(cands, bb, c["best_anchor_off"])
                unpaired_inside = int((c["inside"] & ~c["paired"]).sum())
                crit = {"C1": unpaired_inside == 0 and c["n_pairs"] >= MIN_PAIRS,
                        "C2": c["max_abs_ms"] is not None and c["max_abs_ms"] <= TOL * 1000,
                        "C3": ru_n <= RUNNER_UP_FRAC * c["n_pairs"],
                        "C4": abs(c["ppm"]) <= MAX_PPM}
                pidx = c["pulse_idx"][c["paired"]] + 1  # 1-based pulse numbers
                fits[ph] = (bb, c)
                rep["phases"][ph] = {"block": bb, "n_msgs": c["n_msgs"], "n_inside": c["n_inside"],
                                     "n_pairs": c["n_pairs"], "unpaired_inside": unpaired_inside,
                                     "max_abs_resid_ms": r1(c["max_abs_ms"]), "rms_resid_ms": r1(c["rms_ms"]),
                                     "slope_ppm": r1(c["ppm"]), "offset_s": round(c["a"], 4),
                                     "runner_up_pairs": ru_n, "runner_up": ru_where,
                                     "pair_counts_by_block": {str(b): (cands[b]["n_pairs"] if cands[b] else None) for b in blocks},
                                     "pulse_numbers_first_last": [int(pidx.min()), int(pidx.max())] if len(pidx) else None,
                                     "pulse_numbers": pidx.tolist(),
                                     "resid_ms": [round(float(x) * 1000, 1) for x in c["resid"][c["paired"]]],
                                     "criteria": crit, "pass": all(crit.values())}
            # C5 / joint fit, decision
            dec = {"phases_written": [], "maps": {}}
            P = rep["phases"]
            if sn == 1:
                # sub-01: encoding map, extrapolated to recognition
                if "encode" in P and P["encode"]["pass"] and P["encode"]["block"] == run:
                    bb, c = fits["encode"]
                    m = msgs[(run, "recog")]
                    tr = sorted(m)
                    pred_all = c["a"] + c["b"] * np.array([m[t] for t in tr], float) / 1000.0
                    cr = fits.get("recog", (None, None))[1]
                    if cr is not None and P["recog"]["block"] == run:
                        # recognition pulses paired by the same search; lag vs the encoding map
                        x = cr["msg_s"][cr["paired"]]
                        lag = pulses[cr["pulse_idx"][cr["paired"]]] - (c["a"] + c["b"] * x)
                    else:
                        lag = np.array([np.inf])
                    s1a = {"recog_pulse_minus_mapped_msg_ms_median": r1(np.median(lag) * 1000),
                           "spread_ms": r1(np.ptp(lag) * 1000), "n": int(len(lag)),
                           "pass": bool(np.ptp(lag) <= TOL) and len(lag) == len(tr),
                           "inside": bool(np.all((pred_all >= 0) & (pred_all < dur)))}
                    rep["S1a"] = s1a
                    s1_checks.append(s1a)
                    # where the response-defined image time (EyeLink recording stop - RT) falls after the message
                    rs = np.array(diag[run]["rec_rec_stop"], float)
                    rtt = {e["trial"]: e["rt"] for e in events if e["block"] == run and e["phase"] == "recog" and e["type"] == "trial"}
                    v = np.array([rs[t - 1] - rtt[t] - m[t] for t in tr if t in rtt and not math.isnan(rtt[t])])
                    s1a["stop_minus_rt_minus_msg_ms"] = {"median": r1(np.median(v)), "min": r1(v.min()), "max": r1(v.max())}
                    if s1a["pass"] and s1a["inside"]:
                        dec["phases_written"] = ["recog"]
                        dec["pulse_rows"] = {cr["trials"][i]: float(pulses[cr["pulse_idx"][i]])
                                             for i in np.flatnonzero(cr["paired"])}
                        dec["maps"] = {"recog": (c["a"], c["b"])}
                        dec["label"] = {"recog": "eyelink_clock_from_encoding"}
            else:
                good = [ph for ph in ("encode", "recog") if ph in P and P[ph]["pass"]]
                if len(good) == 2 and P["encode"]["block"] != P["recog"]["block"]:
                    good = []  # phases pointing at different blocks: ambiguous
                    rep["C5"] = {"pass": False, "why": "phases fit different EyeLink blocks"}
                if len(good) == 2:
                    ce, cr = fits["encode"][1], fits["recog"][1]
                    # offset agreement at the boundary between phases (mid-point between last enc, first rec msg)
                    tb = (ce["msg_s"].max() + cr["msg_s"].min()) / 2
                    dd = (ce["a"] + ce["b"] * tb) - (cr["a"] + cr["b"] * tb)
                    c5 = abs(dd) <= TOL
                    x = np.r_[ce["msg_s"][ce["paired"]], cr["msg_s"][cr["paired"]]]
                    y = np.r_[pulses[ce["pulse_idx"][ce["paired"]]], pulses[cr["pulse_idx"][cr["paired"]]]]
                    bj, aj = np.polyfit(x, y, 1)
                    rj = y - (aj + bj * x)
                    rep["C5"] = {"enc_minus_rec_offset_ms": r1(dd * 1000), "pass": bool(c5),
                                 "joint_max_abs_resid_ms": r1(np.max(np.abs(rj)) * 1000),
                                 "joint_rms_ms": r1(np.sqrt(np.mean(rj ** 2)) * 1000), "joint_slope_ppm": r1((bj - 1) * 1e6)}
                    if c5 and np.max(np.abs(rj)) <= TOL:
                        dec["maps"] = {"encode": (aj, bj), "recog": (aj, bj)}
                    else:
                        good = []
                elif len(good) == 1:
                    c = fits[good[0]][1]
                    dec["maps"] = {good[0]: (c["a"], c["b"])}
                dec["phases_written"] = good
                dec["label"] = {ph: "pulse_search_clock_fit" for ph in good}
                if good:
                    dec["block"] = P[good[0]]["block"]
                # S1b: encoding-only extrapolation test on calibration runs of the other participants
                if (sn, run) in MAIN_RUNS and "encode" in P and "recog" in P and P["encode"]["pass"] and P["recog"]["pass"]:
                    ce, cr = fits["encode"][1], fits["recog"][1]
                    pr = ce["a"] + ce["b"] * cr["msg_s"][cr["paired"]]
                    err = pulses[cr["pulse_idx"][cr["paired"]]] - pr
                    receipt["s1b"].append({"participant": f"sub-{sn:02d}", "run": run,
                                           "max_abs_err_ms": r1(np.max(np.abs(err)) * 1000),
                                           "median_err_ms": r1(np.median(err) * 1000)})
            rep["decision"] = {"phases_written": dec["phases_written"], "block": dec.get("block", run if sn == 1 else None),
                               "maps": {k: {"a_s": round(v[0], 5), "slope_ppm": r1((v[1] - 1) * 1e6)} for k, v in dec["maps"].items()}}
            # calibration: authors' indices recovered? agreement with the main table
            if (sn, run) in MAIN_RUNS:
                exp = {"encode": list(range(3, 50, 2)), "recog": list(range(54, 102))}
                rep["recovers_authors_indices"] = {ph: (P[ph]["block"] == run and
                                                        set(P[ph]["pulse_numbers"]) <= set(exp[ph]) and
                                                        P[ph]["n_pairs"] == len(msgs[(run, ph)]))
                                                   for ph in MAIN_RUNS[(sn, run)] if ph in P}
            # rows
            if (sn, run) in TARGETS and dec["phases_written"]:
                blk = dec.get("block", run)
                rows, dropped = build_rows(events, blk, set(dec["phases_written"]), dec["maps"], fs, dur, dec["label"],
                                           dec.get("pulse_rows"))
                outdir = os.path.join(a.out, f"sub-{sn:02d}", "ieeg")
                os.makedirs(outdir, exist_ok=True)
                p = os.path.join(outdir, os.path.basename(vhdr).replace("_ieeg.vhdr", "_events.tsv"))
                write_tsv(p, rows)
                rep["written"] = os.path.relpath(p, a.out)
                rep["n_rows"] = {t: sum(1 for r in rows if r["trial_type"] == t) for t in sorted({r["trial_type"] for r in rows})}
                rep["trials_outside_recording"] = dropped
                rep["onset_range_s"] = [rows[0]["onset"], max(r["onset"] for r in rows)]
            if (sn, run) in MAIN_RUNS and a.main_root and (dec["phases_written"] or sn == 1):
                # agreement with the main table for the phases it holds
                blk = run
                phs = set(dec["phases_written"]) & set(MAIN_RUNS[(sn, run)])
                maps = dict(dec["maps"])
                if sn == 1 and "encode" in fits and P["encode"]["pass"]:
                    maps["encode"] = (fits["encode"][1]["a"], fits["encode"][1]["b"])
                    phs = {"encode"}
                rows, _ = build_rows(events, blk, phs, maps, fs, dur, {p_: "check" for p_ in phs})
                mp = os.path.join(a.main_root, f"sub-{sn:02d}", "ieeg", os.path.basename(vhdr).replace("_ieeg.vhdr", "_events.tsv"))
                if os.path.exists(mp) and phs:
                    L = [l.rstrip("\n").split("\t") for l in open(mp)]
                    h = L[0]
                    main = {(r[h.index("trial_type")], r[h.index("trial")], r[h.index("eyelink_time")]): float(r[h.index("onset")]) for r in L[1:]}
                    ours = {(r["trial_type"], str(r["trial"]), str(r["eyelink_time"])): r["onset"] for r in rows}
                    common = set(main) & set(ours)
                    dd = np.array([ours[k] - main[k] for k in common])
                    rep["agreement_with_main"] = {"rows_main": len(main), "rows_ours": len(ours), "common": len(common),
                                                  "max_abs_diff_ms": r1(np.max(np.abs(dd)) * 1000) if len(dd) else None}
            if dec["phases_written"] and dec["maps"]:
                k0 = dec["phases_written"][0]
                seen = receipt.setdefault("_maps", {}).setdefault(f"sub-{sn:02d}_block-{rep['decision']['block']}", [])
                for (orun, oa, ob, odc, ofs) in seen:
                    # same EyeLink block in two runs: compare the DC channel sample by sample over the overlap
                    a1, b1 = dec["maps"][k0]
                    shift = (oa - a1 * ob / b1)  # other-run time = shift + (ob/b1) * this-run time
                    lag = int(round(shift * fs))
                    best = None
                    for dl in range(-5, 6):
                        i0 = max(0, -(lag + dl)); i1 = min(len(dc), len(odc) - (lag + dl))
                        if i1 - i0 > fs:
                            eq = float(np.mean(dc[i0:i1] == odc[i0 + lag + dl:i1 + lag + dl]))
                            if best is None or eq > best[1]:
                                best = (lag + dl, eq, int(i1 - i0))
                    rep.setdefault("same_block_as", []).append({"run": orun, "other_minus_this_s": round(shift, 3),
                                                                "dc_equal_fraction": best and round(best[1], 6),
                                                                "dc_lag_samples": best and best[0],
                                                                "overlap_samples": best and best[2]})
                seen.append((run, dec["maps"][k0][0], dec["maps"][k0][1], dc, fs))
            receipt["runs"][f"sub-{sn:02d}_run-{run:02d}"] = rep
            print(f"sub-{sn:02d} run-{run:02d}", rep["role"], {ph: (v["block"], v["n_pairs"], v["n_msgs"], v["max_abs_resid_ms"], v["runner_up_pairs"], v["pass"]) for ph, v in P.items()},
                  "->", rep["decision"]["phases_written"], rep.get("S1a"), rep.get("C5"), rep.get("agreement_with_main"), rep.get("recovers_authors_indices"), flush=True)
        # S1 image-onset check: EyeLink stop - rt - message (recognition), per participant
        st = []
        for b in blocks:
            rs = np.array(diag[b]["rec_rec_stop"], float)
            m = msgs[(b, "recog")]
            rt = {e["trial"]: e["rt"] for e in events if e["block"] == b and e["phase"] == "recog" and e["type"] == "trial"}
            st += [rs[t - 1] - rt[t] - m[t] for t in m if t in rt and not math.isnan(rt[t])]
        receipt.setdefault("recog_stop_minus_rt_minus_msg_ms", {})[subj] = {
            "median": r1(np.median(st)) if st else None, "p5": r1(np.percentile(st, 5)) if st else None,
            "p95": r1(np.percentile(st, 95)) if st else None, "n": len(st)}
    s1b = receipt["s1b"]
    receipt["S1b"] = {"n_runs": len(s1b), "max_abs_err_ms": max(x["max_abs_err_ms"] for x in s1b) if s1b else None,
                      "pass": bool(s1b) and max(x["max_abs_err_ms"] for x in s1b) <= TOL * 1000}
    receipt.pop("_maps", None)
    with open(a.receipt, "w") as f:
        json.dump(receipt, f, indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--edfapi", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--receipt", required=True)
    ap.add_argument("--main-root")
    ap.add_argument("--dc-cache")
    main(ap.parse_args())
