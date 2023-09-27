
function parsave_files(str_savefile, variable_type, varargin)

if (variable_type == 1)
    paramsetID = varargin{1}; 
    PP = varargin{2}; 
    save(str_savefile,  'paramsetID', 'PP');

elseif (variable_type == 2)
    psychD = varargin{1}; 
    LG_psychD = varargin{2}; 
    condEpiRes = varargin{3}; 

    save(str_savefile,  'psychD','LG_psychD','condEpiRes' )
    
elseif (variable_type == 3)
        simperf_sub =  varargin{1}; 
        save(str_savefile,  'simperf_sub' )
elseif (variable_type == 4)
        simperf_sub =  varargin{1};
        sim_crit =  varargin{2};

        save(str_savefile,  'simperf_sub', 'sim_crit' )
end


end
