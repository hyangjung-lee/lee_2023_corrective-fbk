
function cube_cond = get_episodecube(clamp_episode)
        for iS = 1 : 5
           tot_trials(iS) = length(clamp_episode{iS}) ; 


           iFb = 2 ; % correct case 
           cube_cond(iS, 1, iFb) = length(find( mod(clamp_episode{iS}, 10) == 0 & sign(clamp_episode{iS} ) == -1)); % correct, d = small 
           cube_cond(iS, 2, iFb) = length(find( mod(clamp_episode{iS}, 10) == 0 & sign(clamp_episode{iS} ) == +1)); % correct, d = large 

            iFb = 1 ; % incorrect case 
           cube_cond(iS, 1, iFb) = length(find( mod(clamp_episode{iS}, 13) == 0 & sign(clamp_episode{iS} ) == -1)); % incorrect, d = small 
           cube_cond(iS, 2, iFb) = length(find( mod(clamp_episode{iS}, 13) == 0 & sign(clamp_episode{iS} ) == +1)); % incorrect, d = large 



        end


end