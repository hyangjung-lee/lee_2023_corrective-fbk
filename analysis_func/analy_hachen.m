%% Episode Analyses on Hachen et al's data

clear all; 

dataind     = 4;  % Hachen et al., 2011; Human data
cd ../

path_home = pwd; 
addpath([path_home,'/lib']); addpath(genpath([path_home,'/analysis_func'])); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_info.histdir = 'prosp'; % run first
input_info.histdir = 'retro'; % run subsequently 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_fold = [path_home, '/data/temp2']; % Results to be saved in this folder 

if strcmp(input_info.histdir , 'prosp')
    str_att = []; 
elseif strcmp(input_info.histdir , 'retro')
    str_att = 'SlowDrift_'; 
end


if dataind == 4
   
     arry = csvread( [path_home,'/data/other_datasets/Hachen_2021_Ncomm/raw_downloaded/Hachen_et_al_2021_Dynamics_of_History_Dependent_Perceptual_Judgment/Human_Dataset.csv'], 1);
     SUB = unique(arry(:,1));
    
    
    for iSS =  1: length(SUB)
        iSub = SUB(iSS); 
        col_trialID = 2; 
        BB = find( arry(:,1) == iSub); 
       
        trial_initialID = BB(  find( arry( BB, col_trialID ) == 1 ) ); 
        numRun = length( trial_initialID ); 
        
       
        kk = 0; 
        for iRun = 1 : numRun
            
           
            if (iRun ~= numRun)
                 vv = [trial_initialID(iRun):     trial_initialID(iRun + 1)- 1 ]' ; 
                 
                 if (length(vv) > 1) % length(vv) = 1 means trial # 1 occurred twice
                     
                     test_presence_dup = find(diff( arry([trial_initialID(iRun):    trial_initialID(iRun + 1)- 1], 2) ) == 0); 
                     
                     if (~isempty(test_presence_dup))
                         
                        dupli_trials = vv( test_presence_dup ) ; 
                        trials_repeat =  [ dupli_trials( find(diff(dupli_trials) ~= 1) ) + 1; dupli_trials(end)]; 
                        lj =[]; last_tidx =  trial_initialID( iRun ) ; 
                        for j = 1 : length(trials_repeat)
                            trials_repset = vv( find( arry([ trial_initialID(iRun):    trial_initialID(iRun + 1)- 1 ], 2) == arry(trials_repeat(j), 2)) ) ; 
                            kj = [ last_tidx:  trials_repset(1) - 1 ]; 
                            lj = [ lj, kj]; 
                            last_tidx = trials_repset(end) ; 

                           if (j == length(trials_repeat))
                               ej = [last_tidx: trial_initialID(iRun + 1)- 1]; 
                               lj = [ lj, ej]; 
                           end
                        end
                     else
                         lj = vv'; 
                         
                     end
                            
                     trial_j = arry(lj, 2); % validate!  
                     trial_st_end = lj; 
                     
           
                 elseif (length(vv) == 1)
                     
                     trial_j = []; % validate!  
                     trial_st_end = []; 
                
                    
                 end
            
        
            else
                
                
                
                vv = [trial_initialID(iRun):     BB(1) + size(arry(find( arry(:,1) == iSub)),1) - 1 ]' ; 
                test_presence_dup = find(diff( arry( [trial_initialID(iRun):     BB(1) + size(arry(find( arry(:,1) == iSub)),1) - 1 ] ,  2) ) == 0) ; 
                  
                 if (~isempty(test_presence_dup)) 
                  
                        dupli_trials = vv( test_presence_dup ); 
                        trials_repeat =  [ dupli_trials( find(diff(dupli_trials) ~= 1) ) + 1; dupli_trials(end)]; 
                        lj =[]; last_tidx = trial_initialID( iRun ) ; 
                        for j = 1 : length(trials_repeat)
                            trials_repset = vv( find( arry([trial_initialID(iRun):   BB(1) + size(arry(find( arry(:,1) == iSub)),1) - 1 ], 2) == arry(trials_repeat(j), 2)) ); 
                           kj = [ last_tidx:  trials_repset(1) - 1 ]; 

                           lj = [ lj, kj]; 

                           last_tidx = trials_repset(end) ; 

                           if (j == length(trials_repeat))
                               ej = [last_tidx: BB(1) + size(arry(find( arry(:,1) == iSub)),1) - 1]; 
                               lj = [ lj, ej]; 
                           end


                        end
                    
                 else
                        lj = vv'; 
                     
                 end
                
                trial_j = arry(lj, 2); % validate!  
                trial_st_end = lj; 
                
                
            end
            col_stimID = 6; 
            stv = unique(  arry(trial_st_end, col_stimID) ); 
            if (isempty(trial_st_end) )
                iRun
            end

            if (~isempty(trial_st_end) )
                kk = kk + 1;  
                iRun
                   Stm{1}{iSS,kk}=  arry(trial_st_end, col_stimID);
                       col_stimID = 5; 

                   Chc{1}{iSS, kk} = arry(trial_st_end, col_stimID); 
                   Chc{1}{iSS, kk}(  Chc{1}{iSS, kk} == 0 ) = -1; 

                   col_stimID = 4; 
                   FBc{1}{iSS, kk} = arry(trial_st_end, col_stimID); 
                   FBc{1}{iSS, kk}(    FBc{1}{iSS, kk}  == 0 ) = - 1; 

            end

            
               
        end
        
        
    end
end

imode = -100
numSim = 1; targetN = 1; flagRT=3; nback = 1;

if (dataind == 3 || dataind == 4)
    clear runSet
    for iCond = 1 
        for iSub = 1 : size(Stm{iCond}, 1)
            k= 0; 
            for iRun = 1 : size(Stm{iCond}, 2)
                if (~ isempty(Stm{iCond}{iSub,iRun}) )

                    k = k+1; 


                end
             end
            runSet(iSub,iCond) = k;  
        end
    end
    Fbvect = 1; 
end
clear MssMap_pre MssMap_post MccMap_pre MccMap_post McfMap_pre McfMap_post
numSub = length(SUB); WhichQuad = SUB; 
nStim = length(unique(Stm{1}{1,1}));
for iG = 1 : 1
    for iSub = 1:length(WhichQuad)

        subIDtoFit = WhichQuad(iSub);
        

         for iFb = Fbvect 
            for iSim = 1 : numSim
                %{
                if imode == 0 
                       ss = subjSmat{iSub, iFb};  ff = subjCLmat{iSub, iFb};
                       cc = subjCmat{iSub, iFb}; 
                end
                if (imode == 100)
                    ss = matT{iSub, iFb}.S;  
                    ff = matT{iSub, iFb}.F;
                    cc = Mwm_DynCrit{iG}.Mchoice{iSub}{iFb}.C(:,:,iSim);
                    finalC =   Mwm_DynCrit{iG}.Mchoice{iSub}{iFb}.critMean(:,:,iSim);
                    temp = zscore(finalC');
                    finalCC = temp'; 
                end
                %}
                if (imode == -100)
                    RR = runSet(iSub, iFb); 
                    for iRun = 1: RR
                        ss = Stm{iFb}{iSub, iRun}'; 
                        CH = Chc{iFb}{iSub, iRun}'; 
                        FB = FBc{iFb}{iSub, iRun}'; 
                        cc = CH; cf = FB; 
                        
pr = 1:  length(ss) - 1; po = 2: length(ss); 
if strcmp(input_info.histdir , 'prosp')
%%%% sorting (Prospective)
    pre = pr;
    post = po;
elseif strcmp(input_info.histdir , 'retro')
     pre = po;
     post = pr;
end
                        for iback = 1 : nback

                            MssMap_pre{iSub}{iback, iFb, iRun}(:, :, iSim) = ss(:, pre(:, 1:end-iback+1));
                            MssMap_post{iSub}{iback, iFb, iRun}(:, :, iSim) = ss(:, post(:, iback:end));      

                            MccMap_pre{iSub}{iback,iFb, iRun}(:,:,iSim)=cc(:,pre(:, 1:end-iback+1));
                            MccMap_post{iSub}{iback,iFb, iRun}(:,:,iSim)=cc(:,post(:, iback:end));

                            McfMap_pre{iSub}{iback, iFb, iRun}(:,:,iSim)=cf(:, pre(:, 1:end-iback+1));
                            McfMap_post{iSub}{iback,iFb, iRun}(:,:,iSim)=cf(:, post(:, iback:end));
%                             if (imode == 100)           
%                                   finalCMap_post{iSub}{iback, iFb, iRun}(:,:,iSim)= finalCC(:, pre(:, 1:end-iback+1));
%                             elseif (imode == 0)
%                                 finalCMap_post{iSub}{iback, iFb, iRun}(:,:,iSim)= ones(size(cf)); 
%                             end
                        end
                    end
                end % imode = -100
              
                
         
            end
         end
    end

    for cond_order = 2 : 7
        psychD.extPure{cond_order} = NaN( numSub, nStim, numSim) ;
        psychD.norPure{cond_order} = NaN( numSub, nStim, numSim) ;
        LG_psychD.extPure{cond_order} = NaN( numSub, nStim, numSim) ;
        LG_psychD.norPure{cond_order} = NaN( numSub, nStim, numSim) ;

        condEpiRes.extPure{cond_order}.SMALL.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.SMALL.numTrials = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numTrials = NaN( numSub, nStim, numSim) ;    
        
        condEpiRes.norPure{cond_order}.SMALL.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.SMALL.numTrials = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numChoices = NaN( numSub, nStim, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numTrials = NaN( numSub, nStim, numSim) ;      
    end
    
    
     
    for iSub = 1 : length(WhichQuad)

        subIDtoFit = WhichQuad(iSub)
        
        
        for iback = targetN:targetN

            ssMap_pre = []; ssMap_post = []; ccMap_pre=[]; ccMap_post=[]; cfMap_pre=[]; cfMap_post = []; crit_pre =[]; crit_post=[];

                        for iFb = Fbvect
                            RR = runSet(iSub, iFb); 
                            DD.ext  = cell(7,  1); LG_DD.ext = cell(7,1);  
                            DD.nor  = cell(7,  1); LG_DD.nor = cell(7,1);  
                        
                            for ii = 2 : 7
                                DD.ext{ii} = cell(1,9);  LG_DD.ext{ii} = cell(1,9); 
                                DD.nor{ii} = cell(1,9);  LG_DD.nor{ii} = cell(1,9); 
                            end
                            
                            for iRun = 1:RR 
                                ssMap_pre= MssMap_pre{iSub}{iback,iFb, iRun}(:,:,iSim )';
                                ssMap_post = MssMap_post{iSub}{iback,iFb, iRun}(:,:,iSim )';

                                ccMap_pre  = MccMap_pre{iSub}{iback,iFb, iRun}(:,:,iSim)';
                                ccMap_post = MccMap_post{iSub}{iback,iFb, iRun}(:,:,iSim)'; 
                                cfMap_pre  = McfMap_pre{iSub}{iback,iFb, iRun}(:,:,iSim)';
                                cfMap_post = McfMap_post{iSub}{iback,iFb, iRun}(:,:,iSim)';


                                                         
                                if (flagRT == 3) 
                                    Xt = ssMap_post;
                                    Xt(Xt <= -4) = -4; 
                                    Xt(Xt > -4 & Xt <= - 3 ) = -3;
                                    Xt(Xt > -3 & Xt <= - 2 ) = -2;
                                    Xt(Xt > -2 & Xt <= - 1 ) = -1;
                                    Xt(Xt >= 1 & Xt <  2 ) = 1;
                                    Xt(Xt >=2  & Xt <  3 ) = 2;
                                    Xt(Xt >= 3 & Xt <  4 ) = 3;
                                    Xt(Xt >= 4  ) = 4;
                                                         

                                    Sprev = ssMap_pre ;   
                                    Chprev = ccMap_pre ; 
                                    Fbprev = cfMap_pre ;                 
                                    if (dataind == 3 || dataind == 4)
                                        range =[-4:4]; S = range;
                                        range = unique(Xt)'; S = range;
                                    else
                                        range = unique(Xt)'; S = range;
                                    end
                                                             
                                     Chpost = ccMap_post ; 
                                     Spost = ssMap_post ; 
                                                       
                                end
                 
                    
                    
               % [2,5,6,7] 2: Weak Stay   /  5: Weak Switch
               % [3,4] non veridical -> now -> 3: Weak Switch / 4: Weak Stay  for EXTRA SIZE
               % XXXS, XXS, XS, S, M , L, XL, XXL, XXXL                                              
                                                           
                                     %%% [1] b2 *
                                     %%% strong-Switch-Condition
                                     %%% -> "Weak Switch"
                                     %%% for EXTRA SIZE
                                        SMALLstrSwM.ext =[];
                                        LARGEstrSwM.ext = []; 
                                        SMALLstrSwM.nor =[];
                                        LARGEstrSwM.nor= []; 


                                        % (1) XXXL, Sch, incorrect -> L more 
                                        idx = find( Sprev == range(9) & Chprev == -1 &  Fbprev == -1); 
                                           
                                        SMALLstrSwM.ext = [SMALLstrSwM.ext; idx];   

                                        % (2) XXL, Sch, incorrect -> L more 
                                        idx = find( Sprev == range(8) & Chprev == -1 &  Fbprev == -1); 
                                        
                                        SMALLstrSwM.nor = [SMALLstrSwM.nor; idx];         



                                        % (3) XXS, Lch, incorrect -> S more 
                                        idx = find( Sprev == range(2) & Chprev == +1 &  Fbprev == -1);
                                          
                                        LARGEstrSwM.nor = [LARGEstrSwM.nor; idx];      

                                        % (4) XXXS, Lch, incorrect -> S more 
                                        idx = find( Sprev == range(1) & Chprev == +1 &  Fbprev == -1); 
                                              
                                        LARGEstrSwM.ext = [LARGEstrSwM.ext; idx];           

                                        %%%%%%%%%%%%%% across all stimuli ###################################
                                        cond_order = 3; 

                                            for iS = 1 : nStim
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.ext); tmp_sensIDX = find(Spost(SMALLstrSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                  tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.nor); tmp_sensIDX = find(Spost(SMALLstrSwM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                            
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.ext); tmp_sensIDX = find(Spost(LARGEstrSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                                
                                                 tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.nor); tmp_sensIDX = find(Spost(LARGEstrSwM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                                                           
                                            end

                                     %%% [2] b3 *
                                     %%% strong-Stay-Condition
                                     %%% -> Weak Stay
                                     %%% EXTRA SIZE
                                       
                                        SMALLstrStyM.ext =[];
                                        LARGEstrStyM.ext =[];

                                        SMALLstrStyM.nor =[];
                                        LARGEstrStyM.nor =[];


                                        % (1) XXXL, Lch, correct -> L more 
                                        idx = find( Sprev == range(9) & Chprev == +1 &  Fbprev == 1);
                                        LARGEstrStyM.ext = [LARGEstrStyM.ext; idx];     


                                        % (2) XXL, Lch, correct -> L more 
                                        idx = find( Sprev == range(8) & Chprev == +1 &  Fbprev == 1);    
                                        LARGEstrStyM.nor = [LARGEstrStyM.nor; idx];    


                                        % (3) XXS, Sch, correct -> L more 
                                        idx = find( Sprev == range(2) & Chprev == -1 &  Fbprev == 1); 
                                        SMALLstrStyM.nor = [SMALLstrStyM.nor; idx];   

                                        % (4) XXXS, Sch, correct -> L more 
                                        idx = find( Sprev == range(1) & Chprev == -1 &  Fbprev == 1); 
                                        SMALLstrStyM.ext = [SMALLstrStyM.ext; idx];   

                                                                
                                        %%%%%%%%%%%%%% across all stimuli ###################################
                                        cond_order = 4; 


                                            for iS = 1 : nStim
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.ext); tmp_sensIDX = find(Spost(SMALLstrStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.nor); tmp_sensIDX = find(Spost(SMALLstrStyM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.ext); tmp_sensIDX = find(Spost(LARGEstrStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.nor); tmp_sensIDX = find(Spost(LARGEstrStyM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                            end
                                            

                                     %%% [3] b4 * weak-Switch-Condition 
                                       
                                        SMALLweakSwM.ext=[];
                                        LARGEweakSwM.ext=[];

                                        SMALLweakSwM.nor=[];
                                        LARGEweakSwM.nor=[];


                                        % (1) XL, Sch, incorrect -> L more 
                                        idx = find( Sprev == range(7) & Chprev == -1 &  Fbprev == -1);                        
                                        SMALLweakSwM.ext = [SMALLweakSwM.ext; idx];      

                                        % (2) L, Sch, incorrect -> L more 
                                        idx = find( Sprev == range(6) & Chprev == -1 &  Fbprev == -1);          
                                        SMALLweakSwM.nor = [SMALLweakSwM.nor;idx];      



                                        % (3) S, Lch, incorrect -> S more 
                                        idx = find( Sprev == range(4) & Chprev == 1 &  Fbprev == -1);
                                        LARGEweakSwM.nor = [LARGEweakSwM.nor; idx]; 

                                        % (4) XS, Lch, incorrect -> S more 
                                        idx = find( Sprev == range(3) & Chprev == 1 &  Fbprev == -1);
                                        LARGEweakSwM.ext = [LARGEweakSwM.ext; idx];


                                        %%%%%%%%%%%%%% across all stimuli ###################################
                                        cond_order = 5; 


                                            for iS = 1 : nStim  
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.ext); tmp_sensIDX = find(Spost(SMALLweakSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.nor); tmp_sensIDX = find(Spost(SMALLweakSwM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.ext); tmp_sensIDX = find(Spost(LARGEweakSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];


                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.nor); tmp_sensIDX = find(Spost(LARGEweakSwM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                        %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                              
                                                
                                            end


                                     %%% [4] b9 * weak-Stay-Condition 
                                      
                                        SMALLweakStyM.ext=[];
                                        LARGEweakStyM.ext=[];
                                        SMALLweakStyM.nor=[];
                                        LARGEweakStyM.nor=[];

                                        % (1) XL, Lch, correct -> L more 
                                        idx = find( Sprev == range(7) & Chprev == 1 &  Fbprev == 1);   
                                        LARGEweakStyM.ext = [LARGEweakStyM.ext; idx]; 


                                        % (2) L, Lch, correct -> L more 
                                        idx = find( Sprev == range(6) & Chprev == 1 &  Fbprev == 1); 
                                        LARGEweakStyM.nor = [LARGEweakStyM.nor; idx];  

                                        % (3) S, Sch, correct -> S more 
                                        idx = find( Sprev == range(4) & Chprev == -1 &  Fbprev == 1);
                                        SMALLweakStyM.nor = [SMALLweakStyM.nor; idx];      

                                        % (4) XS, Sch, correct -> S more 
                                        idx = find( Sprev == range(3) & Chprev == -1 &  Fbprev == 1);
                                        SMALLweakStyM.ext = [SMALLweakStyM.ext; idx];
                                        
                                        %%%%%%%%%%%%%% across all stimuli ###################################
                                        cond_order = 2; 


                                            for iS = 1 : nStim
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.ext); tmp_sensIDX = find(Spost(SMALLweakStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.nor); tmp_sensIDX = find(Spost(SMALLweakStyM.nor) ==  S(iS)); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.ext); tmp_sensIDX = find(Spost(LARGEweakStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.nor); tmp_sensIDX = find(Spost(LARGEweakStyM.nor) ==  S(iS)); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; tmp(tmp_sensIDX)];
                                            %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                            end
                                         
                                    %%% [9] b6 * medium-Swith-Condition 
                           
                                        SMALLmediSwM.ext=[];
                                        LARGEmediSwM.ext=[];


                                        % (1) M, Sch, incorrect -> L more 
                                        idx = find( Sprev == range(5) & Chprev == -1 &  Fbprev == -1);
                                        SMALLmediSwM.ext = [SMALLmediSwM.ext; idx];    

                                        % (2) M, Lch, incorrect -> Smore 
                                        idx = find( Sprev == range(5) & Chprev == 1 &  Fbprev == -1);
                                        LARGEmediSwM.ext = [LARGEmediSwM.ext; idx];    


                                        %%%%%%%%%%%%%% across all stimuli ###################################
                                        cond_order = 7; 

                                            for iS = 1 : nStim
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediSwM.ext); tmp_sensIDX = find(Spost(SMALLmediSwM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                              
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediSwM.ext); tmp_sensIDX = find(Spost(LARGEmediSwM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                         %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            end

                                     %%% [6] b7 * medium-Stay-Condition 

                                           mediStyM = []; 
                                           SMALLmediStyM.ext=[];  LARGEmediStyM.ext=[];
                                           LARGEmediStyM.extL = [];     LARGEmediStyM.extS = []; 
                                           SMALLmediStyM.extL = [];      SMALLmediStyM.extS = []; 

                                             % (1) M, Lch, correct -> L more 

                                             if (dataind == 3 || dataind == 4)
                                                  idx = find( Sprev == range(5) & Chprev == 1 &  Fbprev == 1);
                                             elseif (dataind == 5) 
                                                 if (iFb == 1)
                                                 idx = find( Sprev == range(6) & Chprev == 1 &  Fbprev == 1);
                                                 elseif (iFb == 3)


                                                 idx = find( Sprev == range(4) & Chprev == 1 &  Fbprev == 1);
                                                 end

                                             end
                                                mediStyM = [mediStyM; idx];     
                                                LARGEmediStyM.ext = [LARGEmediStyM.ext; idx];    


                                              % (2) M, Sch, correct -> Smore 

                                                if (dataind == 3 || dataind == 4)
                                                     idx = find( Sprev == range(5) & Chprev == -1 &  Fbprev == 1);
                                                elseif (dataind == 5)
                                                       if (iFb == 1)
                                                             idx = find( Sprev == range(6) & Chprev == -1 &  Fbprev == 1);
                                                       elseif (iFb == 3)
                                                              idx = find( Sprev == range(4) & Chprev == -1 &  Fbprev == 1);
                                                       end
                                                end
                                                mediStyM = [mediStyM; idx];    
                                                SMALLmediStyM.ext = [SMALLmediStyM.ext; idx];    
                                                                
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 6; 


                                            for iS = 1 : nStim   
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediStyM.ext); tmp_sensIDX = find(Spost(SMALLmediStyM.ext) ==  S(iS)); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                           
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediStyM.ext); tmp_sensIDX = find(Spost(LARGEmediStyM.ext) ==  S(iS)); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; tmp(tmp_sensIDX)];
                                                
                                            %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            end
  


                            end

                            
                            
                                      
                                        for cond_order = 2 : 7
                                            for iS = 1 : nStim
                                                psychD.extPure{iG, cond_order}(iSub, iS,  iFb) = sum(DD.ext{cond_order}{iS} == 1)/length(DD.ext{cond_order}{iS}); % DD: episode variable

                                                condEpiRes.extPure{iG, cond_order}.SMALL.numTrials(iSub, iS,  iFb) = length(DD.ext{cond_order}{iS} );
                                                condEpiRes.extPure{iG, cond_order}.SMALL.numChoices(iSub, iS,  iFb) = sum(DD.ext{cond_order}{iS}  == 1);

                                                LG_psychD.extPure{iG, cond_order}(iSub, iS,  iFb) = sum(LG_DD.ext{cond_order}{iS} == 1)/length(LG_DD.ext{cond_order}{iS}); % DD: episode variable
                                                condEpiRes.extPure{iG, cond_order}.LARGE.numTrials(iSub, iS, iFb) = length(LG_DD.ext{cond_order}{iS});
                                                condEpiRes.extPure{iG, cond_order}.LARGE.numChoices(iSub, iS,  iFb) = sum(LG_DD.ext{cond_order}{iS} == 1);

                                                psychD.norPure{iG, cond_order}(iSub, iS,  iFb) = sum(DD.nor{cond_order}{iS} == 1)/length(DD.nor{cond_order}{iS}); % DD: episode variable
                                                condEpiRes.norPure{iG, cond_order}.SMALL.numTrials(iSub, iS,  iFb) = length(DD.nor{cond_order}{iS});
                                                condEpiRes.norPure{iG, cond_order}.SMALL.numChoices(iSub, iS, iFb) = sum(DD.nor{cond_order}{iS} == 1);

                                                LG_psychD.norPure{iG, cond_order}(iSub, iS,  iFb) = sum(LG_DD.nor{cond_order}{iS} == 1)/length(LG_DD.nor{cond_order}{iS}); % DD: episode variable
                                                condEpiRes.norPure{iG, cond_order}.LARGE.numTrials(iSub, iS, iFb) = length(LG_DD.nor{cond_order}{iS});
                                                condEpiRes.norPure{iG, cond_order}.LARGE.numChoices(iSub, iS, iFb) = sum(LG_DD.nor{cond_order}{iS} == 1);



                                            end
                                        end

                        end



             end


                   

                 
                
        

    end

end



mod_str = 'hachen'
iPpair = 0; paramsetID = 0; 
       
mkdir(save_fold)
mkdir( [save_fold, '/', 'obsPsych/' ] ) 

mkdir( [save_fold, '/', 'working/' ] ) 
mkdir( [save_fold, '/', 'subtracted/' ] ) 
mkdir( [save_fold, '/', 'figures/' ] ) 
mkdir( [save_fold, '/', 'backward_pse/' ] ) 
mkdir( [save_fold, '/', 'forward_pse/' ] ) 

input_info.WhichQuad = SUB; 
input_info.dataind = dataind; 

input_info.save_fold = save_fold;
input_info.numSim = numSim;         

input_info.mod_str = mod_str; 
input_info.iPpair = iPpair; 
input_info.paramsetID = paramsetID;

input_info.S = S;
  
if ( isempty( str_att ))   
        epi_file_to_save = strcat( save_fold,'/', 'obsPsych/', 'sim_stats_Ingred', mod_str, '_setID', num2str(iPpair), '_', num2str(numSim),'.mat');
variable_type = 2; 
        parsave_files(epi_file_to_save, variable_type, psychD, LG_psychD, condEpiRes)
        
        input_info.str_att = str_att; 
        input_info.epi_file_to_save = epi_file_to_save; 
       
        OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures( psychD, LG_psychD, condEpiRes, input_info)
end

if (strcmp( str_att, 'SlowDrift_'))   
        epi_file_to_save2 = strcat( save_fold,'/',  'obsPsych/', 'SlowDrift_sim_stats_Ingred', mod_str,  '_setID', num2str(iPpair), '_', num2str(numSim),'.mat');
variable_type = 2; 
        parsave_files(epi_file_to_save2,  variable_type, psychD, LG_psychD, condEpiRes)

        epi_file_to_save = epi_file_to_save2;
        
        input_info.str_att = str_att; 
        input_info.epi_file_to_save = epi_file_to_save;
          
        OtherDataset_ver_par_clean_Xindiv_PSE_saveFigures( psychD, LG_psychD, condEpiRes, input_info)

end
