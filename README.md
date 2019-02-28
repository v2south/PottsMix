# PottsMix

This repository contains software and data associated with the Manuscript “A Potts-Mixture Spatiotemporal Joint Model for Combined MEG and EEG Data". The manuscript could be downloaded from [arXiv](https://arxiv.org/abs/1710.08269) now.

## CONTENTS: 
1.	R_file – R code directory, **PottsMix.R** for running the Iterative Conditional Modes algorithm developed in this paper.
2.	data – directory contains data used in both simulation study and real data application.
3.	result – an empty directory where results will be saved.

**data**  directory contains following R data files and Matlab files:
1.	vert.mat –  the matlab data file for 3D locations of the vertices on cortex. 
2.	Data1 - a subdirectory that contains simulated data ‘Data1.RData’ when K = 3 with two Gaussian signals.
3.	Data2 - a subdirectory that contains simulated data ‘Data2.RData’ when K = 4 with three Gaussian signals.
4.	Data3 - a subdirectory that contains simulated data ‘Data3.RData’ when K = 4 with one sinusoid signal and two Gaussian signals.
5.	Data4 - a subdirectory that contains real data ‘Data4.Data’ used for application. 

6.	For each of ‘Data1.RData’, ‘Data2.RData’ and ‘Data3.RData’, it includes the following data inside:
  * 6.1	X_E - lead field matrix for EEG
  * 6.2	X_M – lead field matrix for MEG
  * 6.3	Y_E – simulated EEG data 
  * 6.4	Y_M – simulated MEG data 
  * 6.5	S – true source of neural activity
  * 6.6	true_Z – true allocation of  mixture components on the cortical surface.
  * 6.7	sim – logical variable indicating if the dataset is simulated(TRUE) or real data(FALSE)

7.	 For ‘Data4.RData’, it only contains following data:
  * 7.1	X_E - lead field matrix for EEG
  * 7.2	X_M – lead field matrix for MEG
  * 7.3	Y_E – real EEG data 
  * 7.4	Y_M – real MEG data 
  * 7.5	sim – logical variable indicating if the dataset is simulated(TRUE) or real data(FALSE)

## Using the software for running ICM algorithm.

1.	Assuming that this folder ‘PottsMix’ has been downloaded to a local disk and unzipped.

2.	Line 12 – 21 contains required R packages used in this algorithm. Please be sure to install all required R packages and their dependencies before running the code. The required R packages are as follows:  *rgl*, *R.matlab*, *scatterplot3d, MASS, PottsUtils, ,MCMCpack, bigRR, nnet, glmnet, fields*.

3.	Line 24 – 37 specifies the dataset that will be used for running this algorithm.  We provide four possibilities. When one particular dataset is used, make sure the other dataset is commented out. To use your own data simply copy it to a data directory with the same format as the datasets we have included.

4.	After you have selected which of the four sample datasets you would like to use in lines 24-37, simply source the entire file and the code will run. After it finishes (a few or so minutes for data of the size included in this package) it will produce graphical output (the plots may take a few minutes to render; the 3D plots can be rotated to see all activation areas), and save results.
If you are running an older version of R you may see some warnings but these can be ignored. Please note that the graphical output will call x11() from the R code. This requires that an X server (such as XQuartz) is installed on your machine. Without an X server installed the code will report an error and the graphics will not display.

