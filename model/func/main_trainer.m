function main_trainer(dataset,imod, mod_str, expInfo, folder_path, mode_lapse)

repWindow = 0.9; 

path_tmp = split(folder_path, '/model'); path_home = path_tmp{1}; 
addpath([path_home, '/analysis_func'])

numSub = expInfo.numSub;
numTrials = expInfo.numTrials;

if (mode_lapse == 1)
    load([path_home, '/data/lapseIndividuals_gl.mat'], 'lapList')
elseif (mode_lapse == 0)
    lapList = zeros(numSub, 1); 
end




% delete(gcp('nocreate')); % parallel processing mode
% parpool('local', numSub);
                
for iSub = 1: numSub
    fitInfo = []; 
    subIDtoFit = iSub; 
    tic;
    fname = [folder_path, '/', 'fitnow_Sub'  num2str(iSub)];
    check_filename = [fname, '.txt']; 
  
    if (exist(check_filename, 'file') == 0)
        fid = fopen(check_filename,'w'); fclose(fid); 
                x0=[];funInputFurther = []; 

               
                funInputFurther.subindx = subIDtoFit;
                funInputFurther.mod_str = mod_str; 
                funInputFurther.imod = imod; 
                funInputFurther.numCond = 3; 
                funInputFurther.numRun = 10; 
                funInputFurther.lapseRate = lapList(iSub); 
                funInputFurther.repWindow = repWindow; 
                funInputFurther.data = dataset; 
            


                tri_id = repmat( [1:numTrials]', funInputFurther.numCond * funInputFurther.numRun, 1);

                des = [];   k = 0;  acrR = []; 
                for iCond = 1 : funInputFurther.numCond
                    for iRun = 1 : funInputFurther.numRun
                        seqRT = funInputFurther.data{subIDtoFit, iCond}.RT(iRun, :);
                        funInputFurther.data{subIDtoFit, iCond}.C(iRun, seqRT > repWindow) = 0; 
                        acrR = [acrR; funInputFurther.data{subIDtoFit, iCond}.C(iRun, 1:numTrials)']; 
                        for iTrial = 1: numTrials
                            k = k+1;
                            des(k,1) = k; 
                            des(k,2) = tri_id(k);
                            des(k,3) = iCond;
                            des(k,4) = iRun;

                        end
                    end
                end

    
    
    
    
    fitInfo.exp_info = des;
    fitInfo.exp_cho = acrR;      
    funInput = funInputFurther; 
    
    % fitting  model to data
    fitInfo.mod_str = mod_str;
    
        if (exist(strcat(folder_path, '/fit_', fitInfo.mod_str , '_subj', num2str(iSub), '.mat'), 'file') == 0)

            if (strcmp(fitInfo.mod_str, 'BMBU'))
                train_bads_BMBU(subIDtoFit, fitInfo, funInput, folder_path)
            elseif (strcmp(fitInfo.mod_str, 'RLVU'))
                train_bads_RLVU(subIDtoFit, fitInfo, funInput, folder_path)
            elseif (strcmp(fitInfo.mod_str, 'HYBR'))
                train_bads_HYBR(subIDtoFit, fitInfo, funInput, folder_path)
            elseif (strcmp(fitInfo.mod_str, 'CFIX'))
                train_bads_FIX(subIDtoFit, fitInfo, funInput, folder_path)
            elseif (strcmp(fitInfo.mod_str, 'ZFIX'))
                train_bads_FIX(subIDtoFit, fitInfo, funInput, folder_path)
            end

        end

    end
end
end

