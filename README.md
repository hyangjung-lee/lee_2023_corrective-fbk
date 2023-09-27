# Lee et al. (2023) PLoS Biology
This repository contains data and processing scripts in support of the preprint:

Lee, H. J., Lee, H., Lim, C. Y., Rhim, I., & Lee, S. H. (2023). What humans learn from corrective feedback for perceptual decision-making. bioRxiv.
https://www.biorxiv.org/content/10.1101/2023.01.11.523567

The paper is currently under revision at PLoS Biology. New link and the doi information will be provided upon publication. 
File names and script description in this repository to indicate figures and tables in the manuscript; they all follow its latest version during the revision process. 
Data and scripts are available at https://osf.io/z57vx/ under a CC-BY 4.0 license.

## Authors
Hyang-Jung Lee, Heeseung Lee, Chae Young Lim, Issac Rhim, and Sang-Hun Lee


## `data`

The [data](data/) folder contains the processed data used to run the analyses and generate the figures.  

The [raw](data/raw/) folder contains raw behavioral data acquired for the current study, which are available both in .csv and [.mat format](data/raw/exp_data.mat). 
- [exp_data.csv](data/raw/exp_data.csv): Data consist of subject ID, condition ID (1: feedback environment manipulated to favor small choice, 2: no favor, 3:  manipulated to favor large choice), run ID, trial ID, stimulus, choice (-1: small choice / 1: large choice), class variable (-1: small class / 1: large class), correct or incorrect feedback (1: correct / 0: incorrect), and response time (sec).  


The [other_datasets](data/other_datasets/) folder contains not only the source data from the following published articles cited in this study but also the processed data for analyses (S7 Fig in Supporting Information). Both datasets were directly downloaded from the respective public repositories:

- [Urai_2017_Ncomm](data/other_datasets/Urai_2017_Ncomm/) folder contains [source data](data/other_datasets/Urai_2017_Ncomm/raw_downloaded), downloaded from [here](http://dx.doi.org/10.6084/m9.figshare.4300043), and the [processed data](data/other_datasets/Urai_2017_Ncomm/Urai_Data_processed/). 

    Urai, A. E., Braun, A., & Donner, T. H. (2017). Pupil-linked arousal is driven by decision uncertainty and alters serial choice bias. Nature communications, 8(1), 14637([link](https://www.nature.com/articles/ncomms14637))


- [Hachen_2021_Ncomm](data/other_datasets/Hachen_2021_Ncomm/) folder contains human [source data](data/other_datasets/Hachen_2021_Ncomm/raw_downloaded), as part of the repository downloaded from [here](https://osf.io/hux4n), and the [processed data](data/other_datasets/Hachen_2021_Ncomm/Hachen_Data_processed/). 

    Hachen, I., Reinartz, S., Brasselet, R., Stroligo, A., & Diamond, M. E. (2021). Dynamics of history-dependent perceptual judgment. Nature communications, 12(1), 6036 ([link](https://www.nature.com/articles/s41467-021-26104-2))

For more information about the structure of the processed data, see [plot_S7Fig.m](plot_S7Fig.m):



The [si_data_plosb](data/si_data_plosb/) folder contains excel spreadsheet files reported as Supporting Information in the paper. All numerical data underlying main and supplementary figures can be found here, along with statistical information mentioned in the Methods section. 
- [S1_Data_lee2023_plosb.xlsx](data/si_data_plosb/S1_Data_lee2023_plosb.xlsx): S1 Data in the paper
- [S2_Data_lee2023_plosb.xlsx](data/si_data_plosb/S2_Data_lee2023_plosb.xlsx): S2 Data in the paper



## `reproducePaper`

To reproduce all figure items that appear in the main text and in Supporting Information, go to the [plot](plot/) folder and run .m files, which also reproduce statistical results reported in the main text as well as in the tables in Supporting Information. 

Each script is executable within the [plot](plot/) folder set as the working directory. 
- plot_Fig#.m - a script to generate one or more figure or table item(s) named after what was referred to in the paper.  

<br/>

## `model`

The [model](model/) folder contains subfolders as follows:

- [main](model/main/): contains scripts to perform model fitting, model simulation, and model recovery. 
To fit the models, go to the [main] folder, run [runme_fit.m](runme_fit.m). To perform model recovery procedure, run [runme_mdrc.m](runme_mdrc.m), which requires model simulations. Model simulator functions (simulator_MODEL.m) can also be found. 

- [func](model/func/): contains scripts needed to train models. The core part of the BMBU model proposed in the paper include [fitting_BMBU.m](model/main/func/fitting_BMBU.m) and [get_BoundLkld.m](model/main/func/get_BoundLkld.m).

- [lib_ver101](model/lib_ver101/): contains functions required to run scripts in the [model](model/) folder 
- bads-master:  [Bayesian Adaptive Direct Search (BADS)](https://github.com/acerbilab/bads) toolbox


Model fitting problems are solved via [Inverse binomial sampling](https://github.com/acerbilab/ibs), in conjunction with BADS. 

- Acerbi, L. & Ma, W. J. (2017). Practical Bayesian Optimization for Model Fitting with Bayesian Adaptive Direct Search. In *Advances in Neural Information Processing Systems 30*, pages 1834-1844. ([link](https://papers.nips.cc/paper/6780-practical-bayesian-optimization-for-model-fitting-with-bayesian-adaptive-direct-search))

- van Opheusden*, B., Acerbi*, L. & Ma, W.J. (2020). Unbiased and efficient log-likelihood estimation with inverse binomial sampling. PLoS Computational Biology 16(12): e1008483. (* equal contribution). ([link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008483))


Bayesian model selection with [Variational Bayes Approach](https://github.com/MBB-team/VBA-toolbox) toolbox (Fig 6D in the paper) 

- Daunizeau, J., Adam, V., & Rigoux, L. (2014). VBA: a probabilistic treatment of nonlinear models for neurobiological and behavioural data. PLoS computational biology, 10(1), e1003441.



<br/>

## `results_exante`

Model behavior to simulate ex ante individual participant's experiments. Parameter sets used to run these simulations can be found in [S1_Data_lee2023_plosb.xlsx](data/si_data_plosb/S1_Data_lee2023_plosb.xlsx). 

## `results_fits`

Model best-fitting fits to individual participant's data

## `results_model_rc`

Model recovery results 

## `results_expost`

Model behavior to simulate ex post individual participant's experiments

## `results_demo`

BMBU's simulated behavior repeating a 6-trial-example-sequence (for demo): relevant to Fig S2 in the paper



<br/>

## `analysis_func`

Codes related to history-dependent (retrospective and prospective) episode analyses in the study:

- [open_episodeFreq.m](analysis_func/open_episodeFreq.m) - script to reproduce results of the episode analyses for data (and for ex post simulations) intended for Fig 7D,E,F
- sort_episodeFreq.m - script to acquire episode info
- [analy_thisset.m](analysis_func/analy_thisset.m) - script to reproduce episode analyses for the current dataset; The [pock](analysis_func/pock) folder has files that are convenient to execute this.

Scripts written to process the data and to utilize functions for loading saved results and running jobs for analyses: 
- [define_lapse_for_individuals.m](analysis_func/define_lapse_for_individuals.m) - a script to define individual lapse rates from grand psychometric curves 
- fit_psych_PSE.m - a script to fit psychometric curves

    Psychometric curve fitting are performed with [Psignifit 4](https://github.com/wichmann-lab/psignifit) located in the psignifit-master folder. 
        Schütt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A. (2016). Painfree and accurate Bayesian estimation of psychometric functions for (potentially) overdispersed data. Vision research, 122, 105-123. ([link](https://www.sciencedirect.com/science/article/pii/S0042698916000390))

- [load_fittedparamInfo.m](analysis_func/load_fittedparamInfo.m) - a function to load individually recovered best-fitting parameters
- [load_gof_models.m](analysis_func/load_gof_models.m) - a function to load metrics of goodness-of-fits 

- get_normCube.m, get_episodecube.m, get_Fpostevent.m - auxilliary functions for the  episode analyses for data


Scripts to reproduce results with the same analyses applied to other datasets:
- [analy_urai.m](analysis_func/analy_urai.m), &nbsp; [analy_hachen.m](analysis_func/analy_hachen.m), &nbsp;   OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures.m,  OtherDataset_plotting_summary_PSE -  auxilliary functions for plotting and saving results





<br/>

## `lib`

Useful open source codes are stored for plotting and statistical tests 

<br/>

## `extra`

Extra information is stored: PSE estimated from an alternative procedure, relevant to Data analysis in the Methods section in the paper. 



## Reference

The paper will be cited here upon publication. 



