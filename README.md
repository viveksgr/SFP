# SFP project scripts

Exhaustive codebase to run the sniffing analyses on the NEMO dataset

## 1. System Requirements

Operating System: Windows

MATLAB Version: R2023B

### Dependencies:
1. [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
2. [GLMSingle](https://github.com/cvnlab/GLMsingle)
3. [Breathmetrics](https://github.com/zelanolab/breathmetrics)

The code was tested on a machine with the following recommended specifications:
RAM: 64GB; processor: AMD Ryzen Threadripper PRO 4.00 GHz or higher;
No non-standard hardware required.

## 2. Installation Guide
### Instructions
To install: 
```bash
git clone https://github.com/viveksgr/SFP.git
```
Install MATLAB Toolboxes directly from the source links provided. 

### Typical Install Time
The installation should take approximately a few minutes on a normal desktop computer.

## 3. Demo
Demonstration scripts to run analyses are provided in the examples folder. The output should be similar to corresponding analyses in Sagar et. al., 2024 "The human brain modulates sniffs according to fine-grained perceptual features of odors" (In preparation)
Runtime for the demo scripts varies by the analyses but most behavioral analyses should be executable within an hour. 

### Instructions to Run on Demo Data
1. Clone the repository and add all dependencies to filepath
2. Launch MATLAB and add <SFP_functions> to path
3. For each demo example, make sure file path for variable <mainroot> is set as that of cloned repository on local machine.
4. Currently, demo files are provided to recreate crucial behavioral results from decoding and multilogistic analyses.

## 4. Instructions for Use and Reproduction:
To repeat the complete analyses using these scripts, acquire the complete raw [dataset](https://www.nature.com/articles/s41593-023-01414-4#data-availability):
The following scripts should be executed in order after preprocessing:
### Behavioral analysis:
1. sfp_behavior_extract: Prepare sniff data across subjects for sniffing analyses
2. sfp_multilogistic: For multilogistic model to predict levels of perceptual descriptor from sniff data
3. sfp_behavioralRSA: For behavioral representational similarity analysis
4. sfp_decoding: For decoding analyses to predict odor identity from sniff shape
5. sfp_timeseries: Time series version of behavioralRSA
### Neural analysis:
1. ARC_createsingle: Create singletrial voxel responses using GLM single package
2. sfp_clustering: To construct ROI based version of clustering RSA
3. sfp_clustering_searchl: To construct searchlight based version of clustering RSA
4. sfp_clustering_searchl_allpercepts_map: To construct searchlight based version of clustering RSA for individual descriptors

*Detailed documentation for all scripts are in progress*

The full analysis may take approximately over 24-48 hours depending on your system specifications and GLM_single settings.

If you encounter issues, contact sgr.vivek@gmail.com.

## 5. Additional Information:
Terms of use: This content is licensed under MIT License.