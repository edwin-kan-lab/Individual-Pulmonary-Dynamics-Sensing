# Pig Study 2024: Respiration 
This repository contains code for the paper titled ['A Near-Field Radio-Frequency System for Continuous Left and Right Lung Volume Sensing' by Kapoor et al](https://ieeexplore.ieee.org/abstract/document/10879308)

## Abstract
Many respiratory disorders are diagnosed by lung function tests, typically performed using spirometers in clinical settings. Non-invasive measures of lung function are of increasing interest to enable convenient, at-home wellness monitoring of respiratory distress. In this work, we used near-field radio-frequency (NFRF) sensors on mechanically ventilated porcine models for continuous acquisition of pulmonary dynamics for each lung independently under various cardiopulmonary interventions. NFRF sensors, previously shown to detect internal tissue motion, were deployed laterally across the pig chest to measure tidal volume during a stepped intervention. A reference spirometer was used for validation. We used bronchial blockers to isolate each lung, allowing us to demonstrate the novel capability of individually monitoring lung volumes. We demonstrated lung volume measurements with an average error of 9.2% between NFRF and reference spirometry. We also showed individual lung dynamics during one-lung obstruction. The results of this work demonstrates NFRF sensors as a novel wearable method for continuous monitoring of individual lung volumes, aiding the management of respiratory distress during acute and chronic conditions, as well as for detecting early onset of life-threatening conditions such as pneumothorax.

## Instructions
The ```code``` folder includes all code files necessary for reproducing all the results reported. Individual interventions can be run by using ```study_main.m``` while all interventions may be run by using ```run_all_exps.m```

The ```manuscript``` folder includes all ```.tex``` files necessary to generate a PDF version of the manuscript as submitted to the journal

## Acknowledgements
This research was supported by the National Institutes of Health under award 1R21EB034562-01.

## Citing Information
If you would like to make use of the code or reference the results in your work, please use the below for citation 
```BibTex
@ARTICLE{10879308,
  author={Kapoor, Aakash and Conroy, Thomas B. and Gangwar, Kapil and Kan, Edwin C. and Araos, Joaquin},
  journal={IEEE Sensors Journal}, 
  title={A Near-Field Radio Frequency System for Continuous Left and Right Lung Volume Sensing}, 
  year={2025},
  volume={25},
  number={7},
  pages={11857-11867},
  keywords={Lungs;Sensors;Monitoring;Wearable devices;Volume measurement;Dielectrics;Belts;Electrical impedance tomography;Antennas;Radio frequency;Lung volumes;near-field radio frequency (NFRF);wearable devices;wellness monitoring},
  doi={10.1109/JSEN.2025.3538501}}
```
