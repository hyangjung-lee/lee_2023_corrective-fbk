clear all; 
%% Results of episode data (script that is executable to reproduce the "Episode Analysis")
cd ../
addpath('./lib'); 
addpath(genpath( './analysis_func')); 
addpath('./model/lib_ver101'); 
path_tmp = split(pwd, '/analysis_func'); path_home = path_tmp{1}; 
path_str.currpath = pwd; 
path_str.data       = [path_home, '/data/raw'];
load([path_str.data '/exp_data.mat'])
nSub = size(subjSmat, 1); 
matT = cell(nSub,3); 
for iSub = 1 : nSub 

    for iCond = 1 : 3
            whichDim =  1; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim); 
            matT{iSub, iCond}.S = Matrix_need;  % stimulus series 10 run x 170 trial

            whichDim =  2; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim); 

            matT{iSub, iCond}.F = Matrix_need;  % feedback series 10 run x 170 trial   

            whichDim =  3; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim);
            matT{iSub, iCond}.C = Matrix_need;  % choice series 10 run x 170 trial

            whichDim =  4; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim);
            matT{iSub, iCond}.RT = Matrix_need;  % choice series 10 run x 170 trial

    end
    
end
clear subjSmat subjCLmat subjCmat subjFmat 
dataset = matT;


numSim = 1; 
numSub = 30; 
modeGroup = 1; 
WhichQuad = 1:numSub; imode = 0; nCond = 3; flagRT = 3; 

input_info.WhichQuad = WhichQuad; 
input_info.modeGroup = modeGroup;
input_info.numSim = numSim; 
input_info.imode = imode;      input_info.numSub = numSub; input_info.nCond = nCond;  input_info.flagRT = flagRT; 


          
          
input_info.histdir = 'retro';
[Pcub] = sort_episodeFreq(matT, subjRTmat, [],input_info); 
toiminusone_Pcub = Pcub; 
input_info.histdir = 'prosp';
[Pcub] = sort_episodeFreq(matT, subjRTmat, [],input_info); 
toiplusone_Pcub = Pcub; 

% save('/Users/eva/Desktop/FF_simul_psinifit/Episode_freq_OBS.mat','toiminusone_Pcub','toiplusone_Pcub'); 
%% same analyses similarly done for ex-post simulation data (scripts written to locally work)
%{
load('Best_param_BMBU.mat')
ParamSetExpost = fit_Param.BMBU; 
ParamSetExpost = [ ParamSetExpost, zeros(30, 1)]; 
simulator_core = @simulator_BMBU; 

imod = 1; 
str_models ={'BMBU', 'RLVU','HYBR'}; 

mod_str = str_models{imod}; 

numSim = 100; 
numSub = 30; 
imode = 100; 
WhichQuad = 1:numSub;
input_info.WhichQuad = WhichQuad; 
input_info.modeGroup = modeGroup;
input_info.numSim = numSim; 
input_info.imode = imode;      input_info.numSub = numSub; input_info.nCond = nCond;
          
          
load('lapseIndividuals_gl.mat', 'lapList'); 
for iSub = 1:numSub
    paramsetID = ParamSetExpost(iSub,:);
    theta_in = paramsetID
    lap_indv = lapList(iSub); 
    matModel = simulate_individual_lapgl(dataset, imod, mod_str, iSub, numSim, theta_in, lap_indv, simulator_core);
    Mwm_DynCrit{1}.Mchoice{iSub} = matModel;
end

input_info.histdir = 'retro';
[Pcub] = sort_episodeFreq(matT, subjRTmat,Mwm_DynCrit,input_info); 
toiminusone_Pcub = Pcub; 

input_info.histdir = 'prosp';
[Pcub] = sort_episodeFreq(matT, subjRTmat, Mwm_DynCrit,input_info); 
toiplusone_Pcub = Pcub; 

save('/Users/eva/Desktop/FF_simul_psinifit/Episode_freq_Wmodel_verS445_FF.mat','toiminusone_Pcub','toiplusone_Pcub'); 


cd '/Users/eva/Desktop/simul_psinifit'


load('Best_param_RLVU.mat')
clear ParamSetExpost
ParamSetExpost = fit_Param.RLVU; 
simulator_core = @simulator_RLVU;



imod = 2; 

mod_str = str_models{imod}; 

 
 for iSub = 1:numSub
    
    paramsetID = ParamSetExpost(iSub,:);
    theta_in = paramsetID;
    lap_indv = lapList(iSub); 
    matModel = simulate_individual_lapgl(dataset, imod, mod_str, iSub, numSim, theta_in, lap_indv, simulator_core);
    Mwm_DynCrit{1}.Mchoice{iSub} = matModel;
end

input_info.histdir = 'retro';
[Pcub] = sort_episodeFreq(matT, subjRTmat,Mwm_DynCrit,input_info); 
toiminusone_Pcub = Pcub; 
input_info.histdir = 'prosp';
[Pcub] = sort_episodeFreq(matT, subjRTmat, Mwm_DynCrit,input_info); 
toiplusone_Pcub = Pcub; 


save('/Users/eva/Desktop/FF_simul_psinifit/Episode_freq_Vmodel_verS_softbias_plus_SSpercpt.mat','toiminusone_Pcub','toiplusone_Pcub'); 

%}
