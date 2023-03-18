# Gene essentiality signature for studying the mechanism of action of drugs

The lack of understanding of a drug's mechanism of action may prevent biomarker identification and ultimately lead to attrition in clinical trials.In this study, we explored whether the integration of loss-of-function genetic and drug sensitivity screening data could define a novel signature to better understand the mechanisms of action of drugs.

# Code and system environment, dependencies and installation guide

Most analysis was conducted using CSC cloud server (cPouta VM) with RstudioServer and R 4.0.5, and should be easy to replicate. To replicate the analysis, load the files from the IO data folder and run the analysis .R or .Rmarkdown files in the corresponding subfolders. The expected output are shown in the .html as rendered notebook output.

The VM used to generate the result has 8 cores and 64GB memory. To replicate the analysis notebook in your own machine, please adjust the input file location and within-machine parallel setting accordingly. 

Heavy parallel computations were conducted using CSC HPC cluster, Puhti, where the R code, bash code and I/O data (described below) is provided. To replicate the analysis in your own cluster, please adapt the bash code for job submission accordingly.

The VM used is operated under Linux distribution CentOS 8.
Note certain within-machine parallel setting are not applicable in Windows system (please adjust multicore parallel mode to multisession parallel if necessary).

The dependencies used in each analysis notebook were recorded by calling sessionInfo() function, result are shown in the rendered Rmarkdown output.To install the needed R packages use install.packages("package_name").  Running time (notebook knitting) was reported in each notebook. 


# Data
The data needed for replicating the analysis are available at 
https://a3s.fi/swift/v1/AUTH_0dc1f1a68ad945fdad8e411c18ef81c1/drug_moa/