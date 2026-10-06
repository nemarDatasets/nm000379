This data set contains intracranial EEG data (ASCII format),
eye-tracking data from an EyeLink 1000 remote system (edf format),
behavioral data, and MATLAB code to reproduce the analyses reported in
the manuscript, “Rapid coordination of effective learning by the human
hippocampus” published in Science Advances.

The file KragelEtal21\_SciAdv.zip contains the raw data
divided into folders according to content type, for each of the six
participants in the study, and the MATLAB code necessary to reproduce
all analyses. MATLAB live scripts provide examples of how to reproduce
the main analyses reported in the manuscript.

External datasets:

In addition to the dataset provided here, three open-access datasets
are analyzed in the manuscript.

	-    The FIGRIM Dataset (http://figrim.mit.edu) 
contains eye-tracking data during a continuous recognition task.

    -    Two additional eye-tracking datasets during free viewing of
repeated scenes are provided in “An
extensive dataset of eye movements during viewing of complex
images,” namely the Memory I and Memory II datasets
(https://datadryad.org/stash/dataset/doi:10.5061/dryad.9pf75).


To reproduce region of interest analyses outside of the hippocampus,
both the seven-network cortical parcellation developed by Yeo,
Krienen et al. (https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011),
and the Harvard-Oxford cortical atlas (https://identifiers.org/neurovault.image:1702)
are required.

Stimuli:

The scenes used in this study are part of Microsoft COCO (https://cocodataset.org).
Scenes were selected from the 2017 Train images. Image identifiers are maintained.

Salience model:

To reproduce analyses that consider the visual salience of each scene,
DeepGaze II model predictions for each stimulus are required. Tensorflow
models and a Jupyter notebook demonstrating their use are available for
download (https://deepgaze.bethgelab.org/).


Software dependencies:

The code in this project was developed using MATLAB r2017b. The
following external packages are required for code execution. Some
external packages are included in the repository.

-   fieldtrip (https://github.com/fieldtrip/fieldtrip)
-   spm12 (https://github.com/spm/spm12)
-   BOSC (https://doi.org/10.1016/j.neuroimage.2010.08.064)
-   Edf2Mat (https://github.com/uzh/edf-converter)
-   boundedline (https://github.com/kakearney/boundedline-pkg)
-   export_fig (https://github.com/altmany/export_fig)

License:

The included code is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or any
later version. See the file COPYING for more details. The release of
this software includes functions from other toolboxes that are covered
under their respective licenses.



