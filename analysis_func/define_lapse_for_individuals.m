
%% define and store individual lapse rates from grand psychometric curves (lapseIndividuals_gl.mat in the /data folder)
clear all; 
path_tmp = split(pwd, '/analysis_func'); path_home = path_tmp{1}; 
path_str.currpath = pwd; 
path_str.data       = [path_home, '/data/raw'];
load([path_str.data '/exp_data.mat'])
addpath(genpath([path_str.currpath, '/psignifit-master']))
nBoot = 5000; 


nT = []; cL = []; 
S = -2:2; 
nCond = size(subjSmat, 2); nSub = size(subjSmat, 1); 
fineTimes = [-2.1:0.001:2.1]';lenAxis = length(fineTimes);
lapsR = zeros(1, nBoot); guessR = zeros(1, nBoot);
for iSub = 1 : nSub
    iSub
  
    for iS = 1 : length(S)
        for iCond = 1:nCond
            nT(:,iS, iCond) = length(find(subjSmat{iSub,iCond} == S(iS) )); 
            cL(:,iS, iCond) = length(find( subjCmat{iSub,iCond}(subjSmat{iSub,iCond} == S(iS) )  == 1) ); % choice large

        end
 
    end
    
       
    nTri = sum(nT, 3);
    nCho = sum(cL, 3);

    yyA = nCho./nTri;
          
          
          
          data(:,1) = S;
          data(:,2) = nCho;
          data(:,3) = nTri; 
          for iboot = 1: nBoot
            [~,stand_aftFit] = fit_psych_PSE(data,  options , iboot) ;
            
              lapsR(iboot) = stand_aftFit(3); 
              guessR(iboot)  = stand_aftFit(4); 
   
          end
  

    lapList(iSub,:) = mean( lapsR ); 
    guessList(iSub,:) = mean( guessR ); 

end

lapList = mean([lapList, guessList],2);

% save([path_home, '/data/lapseIndividuals_gl.mat'],'lapList');
