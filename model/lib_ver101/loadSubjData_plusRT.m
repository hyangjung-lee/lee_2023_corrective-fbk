

function [Matrix_need] = loadSubjData_plusRT( subIDtoFit, iCond, whichDim )

path_tmp = split(pwd, '/model'); path_home = path_tmp{1}; 
data_path = [path_home, '/data/raw/exp_data.mat'];

load(data_path) 
    if (whichDim == 1 )
        Matrix_need = subjSmat{subIDtoFit,iCond}; % stimulus seq.
    end
    if (whichDim == 2 )
          Matrix_need = subjCLmat{subIDtoFit,iCond}; % class seq.
    end
    if (whichDim == 3 )
          Matrix_need = subjCmat{subIDtoFit,iCond}; % choice seq.
    end
     if (whichDim == 4 )
          Matrix_need = subjRTmat{subIDtoFit,iCond}; % rt seq.
    end
end