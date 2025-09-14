# MK-801


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/ScreenNeuroPharm/MK-801/blob/master/LICENSE)

> The repository contains the data and the functions needed to reproduce the analysis reported in the article "In vitro clustered cortical networks reveal NMDA-dependent modulation of repetitive activation sequences".

## Details
All uploaded scripts work with a .mat format. 
To reproduce our analysis is necessary to convert the ```.txt``` format file in ```.mat``` format file using the function ```TxT2Mat.m``` in the folder Conversion. 
All electrophysiological recordings are sampled at 10 KHz. 
```TxT2Mat.m``` function allows obtaining for each electrode (120) the peak train .mat file. 
Peak_train file is a sparse vector that reports the spike occurring, saving the spike amplitude.

### Code folder architecture:

- Conversion folder:
    * Txt2Mat: function to convert ```.txt``` format file in ```.mat``` format file

- DoseResponseCurve folder: 
containes the function to estract the dose response curves


- SpikeAnalysis folder:
    * MFR: function to compute the Mean Firing Rate


- Utilities folder: supplementary functions
