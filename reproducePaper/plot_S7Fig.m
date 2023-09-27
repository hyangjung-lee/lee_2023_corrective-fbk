%% S7 Fig: Run each block to generate each figure panel
% Fig S7A and S7B (based on processed data)
% processing codes can be found; analy_urai.m and analy_hachen.m in the analysis_func folder
%% 1. S7A Fig (from Urai et al.,2021)
clear all; 
cd ../
path_home = pwd; 
addpath(genpath([path_home, '/analysis_func']));
dataind = 11; 

%%%%% toggle conditions : retro, prosp, subtr %%%%%
str_att = []; % 'prosp'
% str_att =  'retro'
% str_att =  'subtract'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([path_home, '/data/other_datasets/Urai_2017_Ncomm/Urai_Data_processed/ExpRDK_Data_Urai.mat']) ;
numSub = size(subjSmat, 1);

S = [20, 35, 45, 48, 52, 55, 65, 80];

mod_str = 'urai'
iPpair = 0; paramsetID = 0;  numSim = 1; 
save_fold = [path_home,  '/data/other_datasets/Urai_2017_Ncomm/Urai_Data_processed']; 

if ( isempty( str_att ))    
    load([save_fold, '/forward_pse/AcSsim_stats_urai_setID0_1.mat'])        
elseif  (strcmp( str_att, 'retro'))
    load([save_fold, '/backward_pse/SlowDrift_AcSsim_stats_urai_setID0_1.mat'])      
elseif  (strcmp( str_att, 'subtract'))    
    load([save_fold, '/subtracted/paramsetID0_1.mat'])      
end
input_info.WhichQuad = 1:numSub; 
input_info.dataind = dataind; 

input_info.save_fold = save_fold;
input_info.numSim = numSim; 
input_info.str_att = str_att; 

input_info.mod_str = mod_str; 
input_info.iPpair = iPpair; 
input_info.paramsetID = paramsetID;

input_info.S = S;
input_info.prefix = str_att;
OtherDataset_plotting_summary_PSE(acSsimIndvsSam, acSsimIndvlSam, acSsimIndvsDiff, acSsimIndvlDiff, input_info)


if ( isempty( str_att ))    
%%% Significance: yellow (toi+1) at stimulus 80% 
    [h, p , ci, stats ] = ttest(acSsimIndvlSam(:,8), 50, 'tail','right');
end

if ( strcmp( str_att, 'subtract' ))    
%%% Significance: blue (subtracted) at stimulus 20% 
    [h, p , ci, stats ] = ttest(acSsimIndvsSam(:,1), 0, 'tail','left');
%%% Significance: yellow (subtracted) at stimulus 80% 
    [h, p , ci, stats ] = ttest(acSsimIndvlSam(:,8), 0, 'tail','right');

end
%% S7B Fig (Human data from Hachen et al.,2021)
% execute, starting from the /reproducePaper folder 
clear all; 
cd ../
path_home = pwd; 
addpath(genpath([path_home, '/analysis_func']));

dataind = 4;
arry = csvread( [path_home, '/data/other_datasets/Hachen_2021_Ncomm/raw_downloaded/Hachen_et_al_2021_Dynamics_of_History_Dependent_Perceptual_Judgment/Human_Dataset.csv'], 1);
column_stim = 6; 
S = unique(arry(~isnan(arry(:,column_stim)),column_stim))'; 

mod_str = 'hachen'
save_fold = [path_home, '/data/other_datasets/Hachen_2021_Ncomm/Hachen_Data_processed']; 
   
   
iPpair = 0; paramsetID = 0; numSim = 1;
SUB = unique(arry(:,1));

clear arry;

%%%%% toggle conditions : retro, prosp, subtr %%%%%
str_att = []; % 'prosp'
% str_att = 'retro'
%  str_att = 'subtract'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( isempty( str_att ))    
    load([save_fold, '/forward_pse/AcSsim_stats_hachen_setID0_1.mat'])        
elseif  (strcmp( str_att, 'retro'))
 
    load([save_fold, '/backward_pse/SlowDrift_AcSsim_stats_hachen_setID0_1.mat'])      
elseif  (strcmp( str_att, 'subtract'))    
    load([save_fold, '/subtracted/paramsetID0_1.mat'])      
end


input_info.WhichQuad = SUB; 

input_info.dataind = dataind; 

input_info.save_fold = save_fold;
input_info.numSim = numSim; 
input_info.str_att = str_att; 

input_info.mod_str = mod_str; 
input_info.iPpair = iPpair; 
input_info.paramsetID = paramsetID;

input_info.S = S;
input_info.prefix = str_att;




OtherDataset_plotting_summary_PSE(acSsimIndvsSam, acSsimIndvlSam, acSsimIndvsDiff, acSsimIndvlDiff, input_info)

if ( isempty( str_att ))    
%%% Significance: blue (toi+1) at stimulus -4 / -3 / -2  
    [h, p , ci, stats ] = ttest(acSsimIndvsSam(:,1), 0, 'tail','left');
    [h, p , ci, stats ] = ttest(acSsimIndvsSam(:,2), 0, 'tail','left');
    [h, p , ci, stats ] = ttest(acSsimIndvsSam(:,3), 0, 'tail','left');
   
end
if ( strcmp( str_att, 'subtract' ))  
    %%% Significance: blue (subtracted) at stimulus -4 

    [h, p , ci, stats ] = ttest(acSsimIndvsSam(:,1), 0, 'tail','left');
end
