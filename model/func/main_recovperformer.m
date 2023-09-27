function main_recovperformer(dataset,imod, mod_str, expInfo, folder_path, mode_lapse)

repWindow = 100; numSim = 1;

path_tmp = split(folder_path, '/model'); path_home = path_tmp{1}; 
addpath([path_home, '/analysis_func'])

numSub = expInfo.numSub;
numTrials = expInfo.numTrials;

if (mode_lapse == 1)
    load([path_home, '/data/lapseIndividuals_gl.mat'], 'lapList')
elseif (mode_lapse == 0)
    lapList = zeros(numSub, 1); 
end


if (imod == 1)
    simulator_core = @simulator_BMBU;
elseif (imod == 2)
    simulator_core = @simulator_RLVU; 
elseif (imod == 3)
    simulator_core = @simulator_HYBR; 
elseif (imod == 4)
    simulator_core = @simulator_fixedNN; 
elseif (imod == 5)
    simulator_core = @simulator_fixedNN;    
end


if (imod == 1)
    fP = load_fittedparamInfo(imod, numSub, path_home);            
    fP(:,2) = 0*ones(size(fP,1), 1);
elseif (imod == 2)
    fP = load_fittedparamInfo(imod, numSub, path_home);            
    fP(:,3) = 0*ones(size(fP,1), 1);
    fP(:,5) = 0.5*ones(size(fP,1), 1);  
elseif (imod == 3)
    fP = load_fittedparamInfo(imod, numSub, path_home);            
    fP(:,2) = 0*ones(size(fP,1), 1);
    fP(:,9) = 0.5*ones(size(fP,1), 1);
elseif (imod == 4)
    fP = load_fittedparamInfo(imod, numSub, path_home);         
elseif (imod == 5)
    fP = load_fittedparamInfo(imod, numSub, path_home);         
    fP(:,2) = 0*ones(size(fP, 1), 1);   
end

% delete(gcp('nocreate')); % parallel processing mode
% parpool('local', numSub);
                
for iSub = 1: numSub
    recovInfo = []; 
    subIDtoFit = iSub; 
    tic;

    if (exist([strcat(folder_path, '/sim_synthetic', mod_str) num2str(iSub), '.mat'], 'file') == 0)
          
                x0=[];funInputFurther = []; 
                funInputFurther.subindx = subIDtoFit;
                funInputFurther.mod_str = mod_str; 
                funInputFurther.imod = imod; 
                funInputFurther.numCond = 3; 
                funInputFurther.numRun = 10; 
                funInputFurther.lapseRate = min(lapList); % same across all subjects 
                funInputFurther.repWindow = repWindow;          


                tri_id = repmat( [1:numTrials]', funInputFurther.numCond * funInputFurther.numRun, 1);
                theta_in = mean( fP , 1); % parameter set in
                
                raw_Simdata  = models_simulate(dataset,  funInputFurther, numSim, theta_in, simulator_core); 
                
                acrR = []; des =[]; 
                for iSim = 1 : numSim
                    k = 0; isimR = []; 
                    for iCond = 1 : funInputFurther.numCond
                        for iRun = 1 : funInputFurther.numRun
                            seqRT = raw_Simdata{iCond}.RT(iRun,:);
                            funInputFurther.data{subIDtoFit, iCond}.S = raw_Simdata{iCond}.S;
                            funInputFurther.data{subIDtoFit, iCond}.F = raw_Simdata{iCond}.F;
                            funInputFurther.data{subIDtoFit, iCond}.C = raw_Simdata{iCond}.C;
                                  funInputFurther.data{subIDtoFit, iCond}.RT = raw_Simdata{iCond}.RT;
                            funInputFurther.data{subIDtoFit, iCond}.C(iRun, seqRT > repWindow, iSim) = 0; 
                            isimR = [isimR; funInputFurther.data{subIDtoFit, iCond}.C(iRun, 1:numTrials, iSim)']; 
                            for iTrial = 1: numTrials
                                k = k+1;
                                des(k,1, iSim) = k; 
                                des(k,2, iSim) = tri_id(k);
                                des(k,3, iSim) = iCond;
                                des(k,4, iSim) = iRun;

                            end
                        end
                    end
                    acrR(:,iSim) = isimR;
                end

        recovInfo.exp_info = des;
        recovInfo.exp_simcho = acrR;                  

        funInput = funInputFurther;                   
        parsave_synthetic([strcat(folder_path, '/sim_synthetic', mod_str) num2str(iSub)], iSub, des, acrR,funInputFurther)

    else
        stored_var = parload_synthetic( [strcat(folder_path, '/sim_synthetic', mod_str) num2str(iSub), '.mat'] ) ; 
        recovInfo.exp_info = stored_var.des;
        recovInfo.exp_simcho = stored_var.acrR;      

        funInput = stored_var.funInputFurther;
    end
    
    
    recovInfo.simulated = mod_str;
 
    

    
    % (1) fitting BMBU model 
    recovInfo.mod_str = 'BMBU';
    if (exist(strcat(folder_path, '/Fit', recovInfo.mod_str , '/Fit', recovInfo.mod_str, '_toSim', recovInfo.simulated, num2str(iSub), '.mat'), 'file') == 0)
         recover_bads_BMBU(subIDtoFit, recovInfo, funInput, folder_path)
    end

    
    % (2) fitting RLVU
     recovInfo.mod_str = 'RLVU';
    if (exist(strcat(folder_path, '/Fit', recovInfo.mod_str , '/Fit', recovInfo.mod_str, '_toSim', recovInfo.simulated, num2str(iSub), '.mat'), 'file') == 0)
        recover_bads_RLVU(subIDtoFit, recovInfo, funInput, folder_path)
    end
     % (3) fitting HYBR
    recovInfo.mod_str = 'HYBR';
    if (exist(strcat(folder_path, '/Fit', recovInfo.mod_str , '/Fit', recovInfo.mod_str, '_toSim', recovInfo.simulated, num2str(iSub), '.mat'), 'file') == 0)
        recover_bads_HYBR(subIDtoFit, recovInfo, funInput, folder_path)
    end
    
    
    % (4) fitting ZFIX
    recovInfo.mod_str = 'ZFIX';
    if (exist(strcat(folder_path, '/Fit', recovInfo.mod_str , '/Fit', recovInfo.mod_str, '_toSim', recovInfo.simulated, num2str(iSub), '.mat'), 'file') == 0)
        recover_bads_FIX(subIDtoFit, recovInfo, funInput, folder_path)
    end
    
    
                 
end
end

