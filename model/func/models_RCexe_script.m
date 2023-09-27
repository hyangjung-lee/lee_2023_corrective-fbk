function models_RCexe_script(imod, mode_lapse, path_dir)


addpath(genpath([path_dir, '/model/bads-master']));

str_models ={'BMBU','RLVU', 'HYBR', 'CFIX', 'ZFIX'}; % ZFIX: model "Base"
mod_str = str_models{imod}; 
bads_folder_path = strcat([path_dir, '/model/out_rc/', 'mainData_num'], mod_str ) ;



mkdir(bads_folder_path)


numSub = 30;
matT = cell(numSub,3); 
for iSub = 1 : numSub 

    for iCond = 1 : 3
           whichDim =  1; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim); 
            matT{iSub, iCond}.S = Matrix_need;  % stimulus series 10 run x 170 trial

            whichDim =  2; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim); 

            matT{iSub, iCond}.F = Matrix_need;  % feedback series (class variable, CL) 10 run x 170 trial   

            whichDim =  3; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim);
            matT{iSub, iCond}.C = Matrix_need;  % choice series 10 run x 170 trial

             whichDim =  4; 
            Matrix_need = loadSubjData_plusRT( iSub, iCond , whichDim);
            matT{iSub, iCond}.RT = Matrix_need;  % choice series 10 run x 170 trial

    end
    
end

numTrials = size(matT{1,  1}.S,2);
expInfo.numSub = numSub;
expInfo.numTrials = numTrials; 





dataset = matT; 
main_recovperformer(dataset, imod, mod_str, expInfo, bads_folder_path, mode_lapse) 
% % delete(gcp('nocreate'))



end
