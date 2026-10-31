# Events for runs the main tables leave out: our alignment, not the authors'

**These event tables use our alignment, not the authors'.** The authors' pipeline (ported in
`code/fix-20261030_events.py`, used for the main `sub-*/ieeg/*_events.tsv` files) assigns sync pulse k to trial k with
fixed indices (study = pulses 3:2:49, test = pulses 54:101, `utils/get_onsets.m`), and that rule fails for 7 runs and
for the test phase of sub-01. For those runs, the depositors of this NEMAR copy (Oct 2026) estimated the mapping between the EyeLink
clock and each run's iEEG clock themselves. The tables are kept here, apart from the raw dataset, and do not replace or
change any main `events.tsv`. The authors' released code does not give a usable alignment for these runs/phases (see
"Relation to the authors' analysis").

## Files
- `sub-*/ieeg/sub-*_task-scenerecognition_run-*_events.tsv`: one table per recovered run, with the same columns and
  event content as the main tables (trials, EyeLink fixations/saccades/blinks, behavioural columns) plus `alignment`.
  Columns: `task-scenerecognition_events.json`. For sub-01 the tables hold the test phase only; its study phase is in
  the main tables.
- `alignment_residuals.tsv`: one row per run and phase (recovered, not recovered, and the 24 calibration runs of the
  main tables): EyeLink block found, message and pair counts, pulse numbers used, residuals, clock slope, runner-up,
  pass/fail per criterion, rows written.
- `alignment_receipt.json`: the full output of the alignment script (per-pair residuals, pulse numbers, pair counts
  per EyeLink block, checks).
- Script: `code/fix-20261031_events_our_alignment.py` in the parent dataset (it imports `code/fix-20261030_events.py`
  for EDF/behaviour parsing and pulse detection). Running it needs the annexed EDF/iEEG content and SR Research's
  EDF Access API (`libedfapi.so`).

## Method
1. Same inputs as the main tables: EyeLink trial messages (`encode`, `recog`) and eye events from the EDF files, the
   behavioural arrays, and the sync pulses of the run's DC channel found as in `get_onsets.m` (`findpeaks`, prominence
   2 µV, one sample before the peak).
2. Pairing by search instead of fixed indices. For each EyeLink block of the participant and each phase, every
   (message i, pulse j) anchor gives a clock offset; each message is paired with the nearest pulse within 50 ms; the
   anchor with most pairs wins (tie: lowest RMS). In the study phase only pre-cue pulses are pairing targets (pulses
   followed by another pulse 0.70-0.85 s later: each study trial has a pre-cue and an image pulse, `get_onsets.m`).
3. A linear clock map `t_run = a + b * t_EyeLink` is fitted by least squares on the pairs (2 refit/re-pair rounds;
   one joint fit when both phases pass and agree).
4. Every event is placed through that map: trial rows at the EyeLink trial message, eye events at their EyeLink start
   time; durations are EyeLink durations. Trials whose mapped message falls outside the recording are not written.
5. sub-01 test phase: its sync pulses do not coincide with the `recog` messages (below), so the map fitted on the same
   run's study phase is applied to the test-phase events.

## Acceptance criteria (stated before the left-out runs were examined)
- C1: every message of the phase whose mapped time lies inside the recording is paired, and at least 10 pairs.
- C2: maximum |residual| <= 50 ms after the fit (the main tables' runs: 6-31 ms).
- C3: unique: the runner-up (any other EyeLink block of the participant in the same phase, or the same block with a
  different pulse assignment, offset more than 100 ms away) pairs at most 50% as many messages as the best.
- C4: clock slope |b - 1| <= 1000 ppm.
- C5: when both phases of a run pass, they point to the same EyeLink block and their offsets agree within 50 ms; one
  joint fit is then used and must also stay within 50 ms.
- sub-01 test phase: C1-C4 on the study phase of the same run, plus S1a (the test-phase pulses, paired by the same
  search, sit at a constant lag from the mapped messages: spread <= 50 ms) and S1b (on the other participants' runs,
  the same study-phase-only map predicts the test-phase pulses within 50 ms).
- Two changes after the first run of the script, before any result for the left-out runs was used: the study-phase
  pre-cue restriction in step 2 (without it, a one-pulse shift onto the image pulses pairs every study message too,
  on every run), and a bug fix in C3 (the runner-up offset was compared with the fitted intercept instead of the best
  anchor offset). The thresholds were not changed.

## Calibration on the 24 runs of the main tables
The same search, run blind on the runs the authors' rule aligns, finds the same EyeLink block and exactly the authors'
pulse indices in all 40 phases (and pairs the sub-01 test-phase messages with pulses 54-101), passes C1-C5 everywhere (max residual 7.1-34.0 ms, slope -35 to
+44 ppm, runner-up at most 12/24 study and 10/48 test pairs) and reproduces the main tables' onsets within 7-34 ms
(every row matched). S1b: across the 16 runs with both phases (sub-02, -04, -05, -06), the study-phase-only map
predicts the test-phase pulses within 32.6 ms.

## Results per run (residuals: max |residual| / RMS, ms; runner-up = pairs of the next-best assignment)
| run | EyeLink block | phase: pairs, residual, runner-up | written |
|---|---|---|---|
| sub-01 run-01..08 | 1..8 | study map: 24/24 pairs, 7.1-13.3 / 3.7-6.1 ms, runner-up 8-11/24; S1a lag 808-1175 ms, spread 15.8-27.6 ms | test phase, 48 trials per run (384) |
| sub-03 run-01 | 1 | study 24/24, 8.8 ms, runner-up 13/24 (block 4): C3 fails; test 48/48, 10.0 / 4.7 ms, runner-up 8/48 | test phase only (48 trials) |
| sub-03 run-02 | 2 | study 24/24, 8.3 ms, runner-up 10/24; test 48/48, 9.4 ms, runner-up 7/48; joint 9.4 / 4.8 ms, phase offsets differ by 1.2 ms | both (24 + 48) |
| sub-03 run-04 | 4 | study 24/24, 11.3 ms, runner-up 13/24 (block 1): C3 fails; test 48/48, 11.3 / 5.5 ms, runner-up 8/48 | test phase only (48) |
| sub-05 run-01 | 1 | study: no pre-cue pulse pairs (the export starts after the study phase); test 48/48 on pulses 3-50, 15.2 / 5.6 ms, runner-up 5/48 | test phase (48) |
| sub-06 run-01 | 2 | study 24/24, 8.4 ms, runner-up 11/24; test 48/48, 14.1 ms, runner-up 7/48; joint 11.5 / 5.3 ms, offsets differ by 4.1 ms | both (24 + 48), block = 2 |
| sub-06 run-02 | 3 | study 24/24, 13.3 ms, runner-up 8/24; test 36/36 inside the recording, 10.0 ms, runner-up 7; joint 12.7 / 5.4 ms | both (24 + 36; test trials 37-48 lie after the end of the file) |
| sub-06 run-04 | 4 | study 14/14 inside the recording (pulses 1-27), 8.0 ms, runner-up 7/14; test 48/48 (pulses 32-79), 23.4 ms, runner-up 5/48; joint 23.8 / 6.1 ms, offsets differ by 3.4 ms | both (14 + 48; study trials 1-10 precede the start of the file) |

Per-phase numbers, clock slopes and the pulse numbers used are in `alignment_residuals.tsv` and
`alignment_receipt.json`.

## Not recovered
- sub-03 runs 1 and 4, study phase: the best mapping pairs all 24 messages (max residual 8.8 and 11.3 ms), but another
  EyeLink block pairs 13 of 24 (54% > 50%), so C3 fails and nothing is written. The sync channel of sub-03 carries
  several peaks per pulse (199-222 instead of 102), which raises chance pairings.
- sub-05 run 1, study phase: not in the recording.
- sub-06 run 2, test trials 37-48 and sub-06 run 4, study trials 1-10: outside the recording.

## Notes
- sub-01 test phase: the `recog` messages precede the paired sync pulse by a block-constant 808-1175 ms, whereas in all
  other participants message and pulse coincide (-4 to +4 ms). The pulse holds the same relation to the trial as in the
  other participants: the response-defined time (EyeLink recording stop minus response time) follows the pulse by
  833-846 ms (run medians), as it follows the message in the other participants (5th-95th percentile 828-847 ms). So
  in sub-01 the `recog/sync_pulse` row (onset of the paired pulse) marks the moment that `recog/trial` marks in every
  other participant, and the `recog/trial` row marks the earlier message. Eye events are placed by the EyeLink clock
  and do not depend on which marker is used.
- sub-06 run 2 is a second export of block 3: its DC channel is sample-identical to run 3's over the whole 398 s
  overlap, shifted by 7.000 s (run 2 starts 7 s before run 3 and ends, 405 s later, inside run 3). Its events are block 3's events on run
  2's time base; do not count them as separate trials.
- sub-06 run 1 holds EyeLink/behavioural block 2 (no iEEG run holds block 1).

## Relation to the authors' analysis
The authors' code does not align these runs: `get_onsets.m` sets sub-06 block 4 onsets to NaN ("missing encoding data
for half of block 4"), uses the fixed indices for the other runs, and applies `adjust_S1_events.m` to sub-01's test
phase, which as released evaluates to 0 for every trial (see the main README). For sub-03, sub-05 run 1 and sub-06
runs 1-2 the fixed indices pair pulses with the wrong trials or block (residuals of seconds or more, listed in the main
README). These tables therefore cannot be checked against the paper's analysis.
