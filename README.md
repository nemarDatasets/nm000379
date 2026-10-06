# Rapid coordination of effective learning by the human hippocampus — intracranial EEG + eye tracking

Stereo-EEG from **6 patients** with medically refractory epilepsy (Northwestern Memorial Hospital Comprehensive Epilepsy
Center, Chicago) performing a scene recognition memory task with concurrent EyeLink eye tracking. Re-packaged in
iEEG-BIDS from the authors' public release (Zenodo record [4728229](https://doi.org/10.5281/zenodo.4728229), v1.0.0,
"Data and scripts related to: Rapid coordination of effective learning by the human hippocampus", CC BY 4.0).

Reference article: Kragel JE, Schuele S, VanHaerents S, Rosenow JM, Voss JL (2021). *Rapid coordination of effective
learning by the human hippocampus.* Science Advances 7(25):eabf7144. https://doi.org/10.1126/sciadv.abf7144

## Participants
Six participants (three females; average age 29 years, range 24–38) with medically refractory epilepsy and depth
electrodes including the hippocampus (inclusion criterion), implanted for neurosurgical monitoring before elective surgery
(Kragel et al. 2021). Sex and age are reported only at group level, so `participants.tsv` lists `n/a` per participant.
The NEMAR dataset on006065 (Kragel et al. 2025, closed-loop theta stimulation) comes from the same laboratory; it is a
different study, task and set of recordings.

## Task
Scene recognition memory task (Presentation 18.0): eight blocks; in each block the participant studied a sequence of 24
images (natural scenes from Microsoft COCO 2017 train images with people, animals or food; eight of each category per
block) followed by a recognition test. Each BIDS run is one block exported by the authors (`S<N>_sl_block<k>.m00` →
`run-<k>`; participant 3 has blocks 1, 2 and 4 only; participants 4–6 fewer blocks, as released). The scene images are
not redistributed (COCO); the image identifiers are in `sourcedata/KragelEtal21_SciAdv/data/behav/S*/S*_stim_array*.txt`.

## Recording
Nihon Kohden amplifier, 1 or 2 kHz per clinical needs (participants 1–3: 2000 Hz; 4–6: 1000 Hz, from the export headers),
hardware band-pass 0.6–600 Hz. Clinical reference and ground: an implanted strip facing the scalp or a scalp electrode.
AD-TECH depth electrodes (contacts 5–10 mm apart). Participant 6 also has scalp channels (`sCz`, `sPz`, typed EEG) and a
`Ref3` channel (MISC).

## What was converted, and how
- Source: Nihon Kohden ASCII exports (`.m00`), header `TimePoints= Channels= BeginSweep[ms]=0.00 SamplingInterval[ms]=…
  Bins/uV=1.000 Time=…`, then one row per sample with microvolt values printed with 2 decimals. Each value was written
  as a BrainVision `INT_32` count with resolution 0.01 µV — **lossless** (the converter asserts that value×100 is an
  integer for every sample, and the round-trip re-reads every ASCII file).
- Channel names: as in the export header with the export's `*` flag characters and a split "`REF 1`" token removed
  (the same clean-up as the authors' `ascii_to_h5.m`); the original header token is kept in `channels.tsv:source_label`.
  The `reference` column is the part after the hyphen (e.g. `REF1`). Participant 3's exports list `J2-REF1` twice
  (columns with different data positions); both are kept, the repeat is named `J2-REF1_2` (noted in
  `status_description`).
- The `DC03`/`DC04` channel is the DC input that the authors' code uses as the **sync-pulse** channel (`ecog_reref.m`);
  typed TRIG. Its header unit label (`(cm)` or empty) is recorded in `status_description`; the values are written in the
  same scaling as all other channels.
- No bipolar re-referencing, line-noise removal or epileptiform-data exclusion (all done in the paper's analysis) was
  applied.
- **Events**: the release has no event table aligned to the iEEG samples. Task events in the paper are built by the
  authors' code from the EyeLink EDF files, the behavioural arrays and the sync pulses (`code/preproc/create_events.m`,
  `adjust_S1_events.m`). These inputs and the code are included under `sourcedata/`, but no `events.tsv` is generated
  here, to avoid inventing an alignment.
- **Electrodes**: contact coordinates from `data/localization/S*/S*_mni.csv` (authors' MNI coordinates: T1 normalized with
  SPM12, deformations applied to CT-identified contacts with Bioimage Suite). They are labelled `space-IXI549Space`, the
  BIDS label for SPM12's normalization template; the extra csv columns are kept (`source_*`).

## sourcedata/
`sourcedata/KragelEtal21_SciAdv/` holds the release content: the original `.m00` files, behavioural arrays, localization
CSVs, MATLAB code (GPL-3.0-or-later per the release; includes third-party toolboxes under their own licences) and the
EyeLink `.edf` files. **EyeLink files were de-identified**: the `** DATE:` header line (recording date and time) was
changed to keep year, month and clock time, with the day set to 01 and the weekday masked (`---`), same byte length;
nothing else was changed (SHA-256 of original and public versions in `sourcedata/b2zen_provenance_IEEG041.json`). The
unmodified originals remain only in the Zenodo release. External datasets used in the paper (FIGRIM, Memory I/II, Yeo
parcellation, Harvard-Oxford atlas, DeepGaze II predictions) are not included; see the release readme.

## Licence
Data: CC BY 4.0 (Zenodo metadata license id `cc-by-4.0`). Code in `sourcedata/`: GNU GPL v3 or later, as stated in the
release.
