function [Pcub] = sort_episodeFreq(matT, subjRTmat, Mwm_DynCrit,input_info)
targetN = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  overall bias measurements & Holding for LARGE 
WhichQuad = input_info.WhichQuad; 
modeGroup = input_info.modeGroup; 
numSim = input_info.numSim; 
imode =  input_info.imode; 
numSub =  input_info.numSub; 
nCond=  input_info.nCond; 
flagRT = input_info.flagRT; 
%%%  overall bias measurements & Holding for LARGE 

nback = 1; % consider 1-back 
S = -2:1:2;

SUB = WhichQuad;

for iG = 1 : modeGroup
    for subIDtoFit = 1:length(WhichQuad)

        iSub = WhichQuad(subIDtoFit);
        
        for iFb = 1: nCond
            for iSim = 1 : numSim

                if imode == 0 
                    ss = matT{iSub, iFb}.S;  
                    ff = matT{iSub, iFb}.F;                  
                    cc = matT{iSub, iFb}.C; 
                end


                if (imode == 100)
                    ss = matT{iSub, iFb}.S;  
                    ff = matT{iSub, iFb}.F;
                    cc = Mwm_DynCrit{iG}.Mchoice{iSub}{iFb}.C(:,:,iSim);
                end

                cf_tmp = (cc == ff); cf = double(cf_tmp);  
                idx = find(cf == 0);
                cf(idx) = -1; 
                
                
if strcmp(input_info.histdir , 'prosp')
%%%% sorting (Prospective)
                pre = 1: length(ss) - 1;
                post = 2: length(ss);
                
elseif strcmp(input_info.histdir , 'retro')
%%%%  sorting (Retrospective) - slow drift   
                pre = 2: length(ss);
                post = 1: length(ss) - 1;
end

                for iback = 1 : nback

                    MssMap_pre{iSub}{iback, iFb}(:, :, iSim) = ss(:, pre(:, 1:end-iback+1));
                    MssMap_post{iSub}{iback, iFb}(:, :, iSim) = ss(:, post(:, iback:end));      

                    MccMap_pre{iSub}{iback,iFb}(:,:,iSim)=cc(:,pre(:, 1:end-iback+1));
                    MccMap_post{iSub}{iback,iFb}(:,:,iSim)=cc(:,post(:, iback:end));

                    McfMap_pre{iSub}{iback, iFb}(:,:,iSim)=cf(:, pre(:, 1:end-iback+1));
                    McfMap_post{iSub}{iback,iFb}(:,:,iSim)=cf(:, post(:, iback:end));
                    
                end
            end
        end

    end

    for cond_order = 2 : 7
        psychD.extPure{cond_order} = NaN( numSub, 5, numSim) ;
        psychD.norPure{cond_order} = NaN( numSub, 5, numSim) ;
        LG_psychD.extPure{cond_order} = NaN( numSub, 5, numSim) ;
        LG_psychD.norPure{cond_order} = NaN( numSub, 5, numSim) ;
        
        condEpiRes.extPure{cond_order}.SMALL.numChoices = NaN( numSub, 5, numSim) ;
        condEpiRes.extPure{cond_order}.SMALL.numTrials = NaN( numSub, 5, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numChoices = NaN( numSub, 5, numSim) ;
        condEpiRes.extPure{cond_order}.LARGE.numTrials = NaN( numSub, 5, numSim) ;
              
        condEpiRes.norPure{cond_order}.SMALL.numChoices = NaN( numSub, 5, numSim) ;
        condEpiRes.norPure{cond_order}.SMALL.numTrials = NaN( numSub, 5, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numChoices = NaN( numSub, 5, numSim) ;
        condEpiRes.norPure{cond_order}.LARGE.numTrials = NaN( numSub, 5, numSim) ;
     
    end
    
    
    numRun = size(ss, 1);     
    for subIDtoFit = 1:length(WhichQuad)

        iSub = WhichQuad(subIDtoFit)
        
        if (ismember(iSub, SUB )) 
            for iback = targetN:targetN

            ssMap_pre = []; ssMap_post = []; ccMap_pre=[]; ccMap_post=[]; cfMap_pre=[]; cfMap_post = []; 

                     for iSim = 1 : numSim                                 
                         DD.ext  = cell(7,  1); LG_DD.ext = cell(7,1); DD.nor  = cell(7,  1); LG_DD.nor = cell(7,1); 
    
                            
                            for ii = 2 : 7
                                
                                DD.ext{ii} = cell(1,5);  LG_DD.ext{ii} = cell(1,5); 
                                DD.nor{ii} = cell(1,5);  LG_DD.nor{ii} = cell(1,5);
 
                            end
                            for iFb = 1 : nCond
                                     for iRun = 1 : numRun 
                                        ssMap_pre= MssMap_pre{iSub}{iback,iFb}(iRun,:,iSim )';
                                        ssMap_post = MssMap_post{iSub}{iback,iFb}(iRun,:,iSim )';

                                        ccMap_pre  = MccMap_pre{iSub}{iback,iFb}(iRun,:,iSim)';
                                        ccMap_post = MccMap_post{iSub}{iback,iFb}(iRun,:,iSim)'; 
                                        cfMap_pre  = McfMap_pre{iSub}{iback,iFb}(iRun,:,iSim)';
                                        cfMap_post = McfMap_post{iSub}{iback,iFb}(iRun,:,iSim)';



                                       if (flagRT == 3) 
                                             Xt = ssMap_post;
                                             Xt(Xt == -2) = 5; Xt(Xt == -1) = 6; Xt(Xt == 0) = 7; Xt(Xt == 1) = 8; Xt(Xt == 2) = 9;             
                                             Xt(Xt == 5) = -1; Xt(Xt == 6) = -.5; Xt(Xt == 7) = 0; Xt(Xt == 8) = .5; Xt(Xt == 9) = 1;

                                             Sprev = ssMap_pre ;
                                             Sprev(Sprev == -2) = 5; Sprev(Sprev == -1) = 6; Sprev(Sprev == 0) = 7; Sprev(Sprev == 1) = 8; Sprev(Sprev == 2) = 9;
                                             Sprev(Sprev == 5) = -1; Sprev(Sprev == 6) = -.5; Sprev(Sprev == 7) = 0; Sprev(Sprev == 8) = .5; Sprev(Sprev == 9) = 1;

                                             Chprev = ccMap_pre ; 
                                             Fbprev = cfMap_pre ; 

                                             range = [-1, -0.5, 0, 0.5, 1]; 
                                             Chpost = ccMap_post ; 
                                             Spost = ssMap_post ; 
                                             Fpost = cfMap_post;

                                       end           

                                            %%% [1] b2 * strong-Switch-Condition 

                                            SMALLstrSwM.ext =[];
                                            LARGEstrSwM.ext = []; 
                                            SMALLstrSwM.nor =[];
                                            LARGEstrSwM.nor= []; 
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEstrSwM.extL = [];     LARGEstrSwM.extS = []; 
                                            LARGEstrSwM.norL = [];     LARGEstrSwM.norS = []; 
                                            SMALLstrSwM.extL = [];      SMALLstrSwM.extS = []; 
                                            SMALLstrSwM.norL = [];      SMALLstrSwM.norS = []; 

                                            % (1) XL, Lch, incorrect -> S more 
                                            idx = find( Sprev == range(5) & Chprev == 1 &  Fbprev == -1); 
 
                                            LARGEstrSwM.ext = [LARGEstrSwM.ext; idx];   
                                            LARGEstrSwM.extL = idx(Chpost(LARGEstrSwM.ext) == 1);
                                            LARGEstrSwM.extS = idx(Chpost(LARGEstrSwM.ext) == -1);
                                            % (2) L, Lch, incorrect -> S more 
                                            idx = find( Sprev == range(4) & Chprev == 1 &  Fbprev == -1); 
 
                                            LARGEstrSwM.nor = [LARGEstrSwM.nor; idx];         
                                            LARGEstrSwM.norL = idx(Chpost(LARGEstrSwM.nor) == 1);
                                            LARGEstrSwM.norS = idx(Chpost(LARGEstrSwM.nor) == -1);


                                            % (3) S, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(2) & Chprev == -1 &  Fbprev == -1);

                                            SMALLstrSwM.nor = [SMALLstrSwM.nor; idx];      
                                            SMALLstrSwM.norL = idx(Chpost(SMALLstrSwM.nor) == 1);
                                            SMALLstrSwM.norS = idx(Chpost(SMALLstrSwM.nor) == -1);

                                            % (4) XS, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(1) & Chprev == -1 &  Fbprev == -1); 
    
                                            SMALLstrSwM.ext = [SMALLstrSwM.ext; idx];           
                                            SMALLstrSwM.extL = idx(Chpost(SMALLstrSwM.ext) == 1);
                                            SMALLstrSwM.extS = idx(Chpost(SMALLstrSwM.ext) == -1);
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 3; 

                                            for iS = 1 : 5
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.ext); tmp_sensIDX = find(Spost(SMALLstrSwM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLstrSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];         
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrSwM.nor); tmp_sensIDX = find(Spost(SMALLstrSwM.nor) ==  S(iS) ); 
                                                tmp2 = []; tmp2 = Fpost(SMALLstrSwM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; vect_full];
                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.ext); tmp_sensIDX = find(Spost(LARGEstrSwM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEstrSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];
                                                                                        
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrSwM.nor); tmp_sensIDX = find(Spost(LARGEstrSwM.nor) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEstrSwM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; vect_full];             
                                            end

                                            %%% [2] b3 * strong-Stay-Condition 

                                            SMALLstrStyM.ext =[];
                                            LARGEstrStyM.ext =[];
                                            SMALLstrStyM.nor =[];
                                            LARGEstrStyM.nor =[];
                                            LARGEstrStyM.extL = [];     LARGEstrStyM.extS = []; 
                                            LARGEstrStyM.norL = [];     LARGEstrStyM.norS = []; 
                                            SMALLstrStyM.extL = [];      SMALLstrStyM.extS = []; 
                                            SMALLstrStyM.norL = [];      SMALLstrStyM.norS = []; 

                                            % (1) XL, Sch, correct -> S more 
                                            idx = find( Sprev == range(5) & Chprev == -1 &  Fbprev == 1);

                                            SMALLstrStyM.ext = [SMALLstrStyM.ext; idx];     
                                            SMALLstrStyM.extL = idx(Chpost(SMALLstrStyM.ext) == 1);
                                            SMALLstrStyM.extS = idx(Chpost(SMALLstrStyM.ext) == -1);

                                            % (2) L, Sch, correct -> S more 
                                            idx = find( Sprev == range(4) & Chprev == -1 &  Fbprev == 1); 
 
                                            SMALLstrStyM.nor = [SMALLstrStyM.nor; idx];    
                                            SMALLstrStyM.norL = idx(Chpost(SMALLstrStyM.nor) == 1);
                                            SMALLstrStyM.norS = idx(Chpost(SMALLstrStyM.nor) == -1);

                                            % (3) S, Lch, correct -> L more 
                                            idx = find( Sprev == range(2) & Chprev == 1 &  Fbprev == 1); 

                                            LARGEstrStyM.nor = [LARGEstrStyM.nor; idx];   
                                            LARGEstrStyM.norL = idx(Chpost(LARGEstrStyM.nor) == 1);
                                            LARGEstrStyM.norS = idx(Chpost(LARGEstrStyM.nor) == -1);

                                            % (4) XS, Lch, correct -> L more 
                                            idx = find( Sprev == range(1) & Chprev == 1 &  Fbprev == 1); 
   
                                            LARGEstrStyM.ext = [LARGEstrStyM.ext; idx];   
                                            LARGEstrStyM.extL = idx(Chpost(LARGEstrStyM.ext) == 1);
                                            LARGEstrStyM.extS = idx(Chpost(LARGEstrStyM.ext) == -1); 
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 4; 


                                            for iS = 1 : 5
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.ext); tmp_sensIDX = find(Spost(SMALLstrStyM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLstrStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 

                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLstrStyM.nor); tmp_sensIDX = find(Spost(SMALLstrStyM.nor) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLstrStyM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; vect_full];



                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.ext); tmp_sensIDX = find(Spost(LARGEstrStyM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEstrStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];


                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEstrStyM.nor); tmp_sensIDX = find(Spost(LARGEstrStyM.nor) ==  S(iS));
                                                tmp2 = []; tmp2 = Fpost(LARGEstrStyM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; vect_full];
  
                                                                                           
                                                   

                                            end
                                            

                                            %%% [3] b4 * weak-Switch-Condition  
                                            SMALLweakSwM.ext=[];
                                            LARGEweakSwM.ext=[];
                                            SMALLweakSwM.nor=[];
                                            LARGEweakSwM.nor=[];
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEweakSwM.extL = [];     LARGEweakSwM.extS = []; 
                                            LARGEweakSwM.norL = [];     LARGEweakSwM.norS = []; 
                                            SMALLweakSwM.extL = [];      SMALLweakSwM.extS = []; 
                                            SMALLweakSwM.norL = [];      SMALLweakSwM.norS = []; 

                                            % (1) XL, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(5) & Chprev == -1 &  Fbprev == -1); 
    
                                            SMALLweakSwM.ext = [SMALLweakSwM.ext; idx];      
                                            SMALLweakSwM.extL = idx(Chpost(SMALLweakSwM.ext) == 1);
                                            SMALLweakSwM.extS = idx(Chpost(SMALLweakSwM.ext) == -1);

                                            % (2) L, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(4) & Chprev == -1 &  Fbprev == -1); 
   
                                            SMALLweakSwM.nor = [SMALLweakSwM.nor;idx];      
                                            SMALLweakSwM.norL = idx(Chpost(SMALLweakSwM.nor) == 1);
                                            SMALLweakSwM.norS = idx(Chpost(SMALLweakSwM.nor) == -1);

                                            % (3) S, Lch, incorrect -> S more 
                                            idx = find( Sprev == range(2) & Chprev == 1 &  Fbprev == -1);

   
                                            LARGEweakSwM.nor = [LARGEweakSwM.nor; idx]; 
                                            LARGEweakSwM.norL = idx(Chpost(LARGEweakSwM.nor) == 1);
                                            LARGEweakSwM.norS = idx(Chpost(LARGEweakSwM.nor) == -1);

                                            % (4) XS, Lch, incorrect -> S more 
                                            idx = find( Sprev == range(1) & Chprev == 1 &  Fbprev == -1);
  
                                            LARGEweakSwM.ext = [LARGEweakSwM.ext; idx];
                                            LARGEweakSwM.extL = idx(Chpost(LARGEweakSwM.ext) == 1);
                                            LARGEweakSwM.extS = idx(Chpost(LARGEweakSwM.ext) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 5; 


                                            for iS = 1 : 5   
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.ext); tmp_sensIDX = find(Spost(SMALLweakSwM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLweakSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 

                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakSwM.nor); tmp_sensIDX = find(Spost(SMALLweakSwM.nor) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLweakSwM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.ext); tmp_sensIDX = find(Spost(LARGEweakSwM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEweakSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2);
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakSwM.nor); tmp_sensIDX = find(Spost(LARGEweakSwM.nor) ==  S(iS));
                                                tmp2 = []; tmp2 = Fpost(LARGEweakSwM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2);
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; vect_full];
           
                                            end


                                            %%% [4] b5 * weak-Stay-Condition 
                                            SMALLweakStyM.ext=[];
                                            LARGEweakStyM.ext=[];
                                            SMALLweakStyM.nor=[];
                                            LARGEweakStyM.nor=[];
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            LARGEweakStyM.extL = [];     LARGEweakStyM.extS = []; 
                                            LARGEweakStyM.norL = [];     LARGEweakStyM.norS = []; 
                                            SMALLweakStyM.extL = [];      SMALLweakStyM.extS = []; 
                                            SMALLweakStyM.norL = [];      SMALLweakStyM.norS = []; 

                                            % (1) XL, Lch, correct -> L more 
                                            idx = find( Sprev == range(5) & Chprev == 1 &  Fbprev == 1); 

                                            LARGEweakStyM.ext = [LARGEweakStyM.ext; idx]; 
                                            LARGEweakStyM.extL = idx(Chpost(LARGEweakStyM.ext) == 1);
                                            LARGEweakStyM.extS = idx(Chpost(LARGEweakStyM.ext) == -1);

                                            % (2) L, Lch, correct -> L more 
                                            idx = find( Sprev == range(4) & Chprev == 1 &  Fbprev == 1); 
    
                                            LARGEweakStyM.nor = [LARGEweakStyM.nor; idx];  
                                            LARGEweakStyM.norL = idx(Chpost(LARGEweakStyM.nor) == 1);
                                            LARGEweakStyM.norS = idx(Chpost(LARGEweakStyM.nor) == -1);

                                            % (3) S, Sch, correct -> S more 
                                            idx = find( Sprev == range(2) & Chprev == -1 &  Fbprev == 1);
  
                                            SMALLweakStyM.nor = [SMALLweakStyM.nor; idx];      
                                            SMALLweakStyM.norL = idx(Chpost(SMALLweakStyM.nor) == 1);
                                            SMALLweakStyM.norS = idx(Chpost(SMALLweakStyM.nor) == -1);

                                            % (4) XS, Sch, correct -> S more 
                                            idx = find( Sprev == range(1) & Chprev == -1 &  Fbprev == 1);
  
                                            SMALLweakStyM.ext = [SMALLweakStyM.ext; idx];
                                            SMALLweakStyM.extL = idx(Chpost(SMALLweakStyM.ext) == 1);
                                            SMALLweakStyM.extS = idx(Chpost(SMALLweakStyM.ext) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 2; 


                                            for iS = 1 : 5  
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.ext); tmp_sensIDX = find(Spost(SMALLweakStyM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLweakStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLweakStyM.nor); tmp_sensIDX = find(Spost(SMALLweakStyM.nor) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLweakStyM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                
                                                DD.nor{cond_order}{iS} = [DD.nor{cond_order}{iS}; vect_full];
           
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.ext); tmp_sensIDX = find(Spost(LARGEweakStyM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEweakStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];
      
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEweakStyM.nor); tmp_sensIDX = find(Spost(LARGEweakStyM.nor) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEweakStyM.nor); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.nor{cond_order}{iS} = [LG_DD.nor{cond_order}{iS}; vect_full];                                                         

                                            end
                                         
                                            %%% [5] b6 * medium-Swith-Condition 

                                            SMALLmediSwM.ext=[]; LARGEmediSwM.ext=[];
                                            LARGEmediSwM.extL = [];     LARGEmediSwM.extS = []; 
                                            SMALLmediSwM.extL = [];      SMALLmediSwM.extS = []; 

                                            % (1) M, Sch, incorrect -> L more 
                                            idx = find( Sprev == range(3) & Chprev == -1 &  Fbprev == -1);

                                            SMALLmediSwM.ext = [SMALLmediSwM.ext; idx];    
                                            SMALLmediSwM.extL = idx(Chpost(SMALLmediSwM.ext) == 1);
                                            SMALLmediSwM.extS = idx(Chpost(SMALLmediSwM.ext) == -1);

                                            % (2) M, Lch, incorrect -> Smore 
                                            idx = find( Sprev == range(3) & Chprev == 1 &  Fbprev == -1);

                                            LARGEmediSwM.ext = [LARGEmediSwM.ext; idx];    
                                            LARGEmediSwM.extL = idx(Chpost(LARGEmediSwM.ext) == 1);
                                            LARGEmediSwM.extS = idx(Chpost(LARGEmediSwM.ext) == -1);

                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 7; 
                                            
                                            for iS = 1 : 5
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediSwM.ext); tmp_sensIDX = find(Spost(SMALLmediSwM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(SMALLmediSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 

                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediSwM.ext); tmp_sensIDX = find(Spost(LARGEmediSwM.ext) ==  S(iS)); 

                                                tmp2 = []; tmp2 = Fpost(LARGEmediSwM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 

                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];    
                                         %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                                            end

                                         %%% [6] b7 * medium-Stay-Condition 

                                            SMALLmediStyM.ext=[];  LARGEmediStyM.ext=[];
                                            LARGEmediStyM.extL = [];     LARGEmediStyM.extS = []; 
                                            SMALLmediStyM.extL = [];      SMALLmediStyM.extS = []; 

                                            % (1) M, Lch, correct -> L more 
                                            idx = find( Sprev == range(3) & Chprev == 1 &  Fbprev == 1);
    
                                            LARGEmediStyM.ext = [LARGEmediStyM.ext; idx];    
                                            LARGEmediStyM.extL = idx(Chpost(LARGEmediStyM.ext) == 1);
                                            LARGEmediStyM.extS = idx(Chpost(LARGEmediStyM.ext) == -1);

                                            % (2) M, Sch, correct -> Smore 
                                            idx = find( Sprev == range(3) & Chprev == -1 &  Fbprev == 1);
  
                                            SMALLmediStyM.ext = [SMALLmediStyM.ext; idx];    
                                            SMALLmediStyM.extL = idx(Chpost(SMALLmediStyM.ext) == 1);
                                            SMALLmediStyM.extS = idx(Chpost(SMALLmediStyM.ext) == -1);
                                            %%%%%%%%%%%%%% across all stimuli ###################################
                                            cond_order = 6; 


                                            for iS = 1 : 5   
                                              

                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(SMALLmediStyM.ext); tmp_sensIDX = find(Spost(SMALLmediStyM.ext) ==  S(iS)); 
                                                
                                                tmp2 = []; tmp2 = Fpost(SMALLmediStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                DD.ext{cond_order}{iS} = [DD.ext{cond_order}{iS}; vect_full];
                                                
                                                tmp = []; tmp_sensIDX =[]; tmp = Chpost(LARGEmediStyM.ext); tmp_sensIDX = find(Spost(LARGEmediStyM.ext) ==  S(iS)); 
                                                tmp2 = []; tmp2 = Fpost(LARGEmediStyM.ext); 
                                                [vect_full] = get_Fpostevent(tmp_sensIDX, tmp, tmp2); 
                                                LG_DD.ext{cond_order}{iS} = [LG_DD.ext{cond_order}{iS}; vect_full];                                      
                                            %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            end
  


                                     end
                            end
                            for cond_order = 2 : 7
                                clamp_episode = DD.ext{cond_order}; cube_cond = []; 
                                cube_cond = get_episodecube(clamp_episode)   ;    
                                cub.extPure{cond_order}.SMALL(:,:,:, iSim) = cube_cond;

                                clamp_episode = DD.nor{cond_order}; cube_cond = []; 
                                cube_cond = get_episodecube(clamp_episode)   ;
                                cub.norPure{cond_order}.SMALL(:,:,:, iSim) = cube_cond;


                                clamp_episode = LG_DD.ext{cond_order}; cube_cond = []; 
                                cube_cond = get_episodecube(clamp_episode)   ;    
                                cub.extPure{cond_order}.LARGE(:,:,:, iSim) = cube_cond;

                                clamp_episode = LG_DD.nor{cond_order}; cube_cond = []; 
                                cube_cond = get_episodecube(clamp_episode)   ;
                                cub.norPure{cond_order}.LARGE(:,:,:, iSim) = cube_cond;
                            end

                     end

             end
        
        end
        
        for cond_order = 2 : 7                                 
            Pcub.extPure{cond_order}.SMALL(:,:,:, iSub) = sum( cub.extPure{cond_order}.SMALL, 4 );
            Pcub.norPure{cond_order}.SMALL(:,:,:, iSub) =  sum( cub.norPure{cond_order}.SMALL, 4 );
            Pcub.extPure{cond_order}.LARGE(:,:,:, iSub) =  sum( cub.extPure{cond_order}.LARGE, 4 );
            Pcub.norPure{cond_order}.LARGE(:,:,:, iSub) = sum( cub.norPure{cond_order}.LARGE, 4 );     
        end                        

    end
     
end

clear DD MssMap_pre MssMap_post MccMap_pre MccMap_post McfMap_pre McfMap_post 

%% Pcub = cub = condEpiRes; Out variable containing "all sorted episodes"
     
end
