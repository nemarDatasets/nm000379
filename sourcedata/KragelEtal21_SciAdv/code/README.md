# Rapid coordination of effective learning by the human hippocampus

This repository is the official implementation of [Rapid coordination of effective learning by the human hippocampus](https://doi.org/10.1101/2020.10.20.347831).

### Dependencies

The code in this projects was developed using MATLAB r2017b. The following external packages are required for code execution. Some external packages are included in the repository.

- [fieldtrip](https://github.com/fieldtrip/fieldtrip)
- [spm12](https://github.com/spm/spm12)
- Better oscillation detection ([BOSC](https://doi.org/10.1016/j.neuroimage.2010.08.064))
- [Edf2Mat](https://github.com/uzh/edf-converter) toolbox for reading EyeLink EDF data.
- [boundedline](https://github.com/kakearney/boundedline-pkg) 
- [export_fig](https://github.com/altmany/export_fig)
- [MES](https://github.com/hhentschke/measures-of-effect-size-toolbox) toolbox for computing effect sizes.

### Datasets

To reproduce all results presented in the manuscript, three datasets are required for download. The primary data is hosted on [Zenodo](10.5281/zenodo.4728229). In addition, three independent datasets are analyzed. The first is the [FIGRIM Dataset ](http://figrim.mit.edu/)which contains [eye-tracking data](http://figrim.mit.edu/Targets.zip) during a continuous recognition task. Two additional eye-tracking datasets during free viewing of repeated scenes are provided in [An extensive dataset of eye movements during viewing of complex images](https://datadryad.org/stash/dataset/doi:10.5061/dryad.9pf75), namely the Memory I and Memory II datasets.

To reproduce region of interest analyses outside the hippocampus, both the seven network cortical parcellation developed by [Yeo, Krienen et al.](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011) (available for download [here](ftp://surfer.nmr.mgh.harvard.edu/pub/data/Yeo_JNeurophysiol11_MNI152.zip)) and the [Harvard-Oxford cortical atlas](https://identifiers.org/neurovault.image:1702).

### Stimuli

The scenes used in this study are part of [Microsoft COCO: Common Objects in Context](https://cocodataset.org). Scenes were selected from the [2017 Train images](http://images.cocodataset.org/zips/train2017.zip). Image identifiers are maintained. 

### Salience model

In order to reproduce analyses that consider the visual salience of each scene, [DeepGaze II](https://deepgaze.bethgelab.org/) model predictions for each stimulus are required. Tensorflow models and a Jupyter notebook demonstrating their use are available for [download](https://drive.google.com/open?id=1kYUwoatqQUS5EabeeSDc6gRmCysnVZ6N).

### Using the software

The code is organized by analysis type. Documentation on how to reproduce the main figures presented in the manuscript is provided within each analysis directory (e.g., [the primary analysis of theta oscillations](link)). In general, a few high-level wrappers will reproduce the main results.

### License

This project is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. See the file COPYING for more details.