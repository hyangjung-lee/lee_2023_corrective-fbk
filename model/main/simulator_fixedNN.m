%% Fixed model simulator (CFIX ZFIX) 
function [matModel]= simulator_fixedNN( theta,~ , InputFurther)
 
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx; 
mod_str = InputFurther.mod_str; 
imod = InputFurther.imod; 

numCond =  InputFurther.numCond; 
numSim = InputFurther.numSim; 
 
numBlock_perRun   = 34;
numRun            = InputFurther.numRun;
numTrial_perBlock = 5;

tic;
      
  
mSig =  theta(1); 
z_mu_init = theta(2);
lapse = InputFurther.lapseRate;



% mu0 = 0; ww = 10;
% ran_mu = 10;resol_mu = 200/ww; % resolution of mu
% z_range = linspace(mu0 - ran_mu, mu0 + ran_mu, resol_mu); len_ax = length(z_range); 
% Prob = nan(1,len_ax);

numTrial_perRun = numBlock_perRun*numTrial_perBlock;
   
for iSim = 1 : numSim    
    for iCond = 1:numCond
        
        seqS = matT{subIDtoFit, iCond}.S(1:numRun, : ); wh_size =  size(seqS); pC_vect = zeros(wh_size);
        seqF = matT{subIDtoFit, iCond}.F(1:numRun, : );     
        lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 
        for iRun = 1 : size(seqS, 1)
           lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
           lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
        end

m_rnd = normrnd(0, mSig, numRun, numTrial_perRun);        
     
            for iTrial = 1: numTrial_perRun

                
                    
x_mu = seqS(:,iTrial);
m = x_mu + m_rnd(:, iTrial);



%% Make choices


cm = z_mu_init*ones(size(x_mu)); 

pC_vect( m - cm > 0 , iTrial) = 1;
pC_vect( m - cm < 0  , iTrial) = -1;
pC_vect( m == cm  , iTrial) = 0 ;

                   

if (sum(lapse_idx(:,iTrial)) ~= 0)
    lap_trial_idx = lapse_idx(:,iTrial)==1; 
    pC_vect(lap_trial_idx, iTrial) = lapse_val(lap_trial_idx,iTrial); 
end


  
            end

           

matModel{iCond}.C(:,:,iSim) = pC_vect;
matModel{iCond}.S = seqS;
matModel{iCond}.F = seqF;
matModel{iCond}.RT = matT{subIDtoFit, iCond}.RT;
                       
            
                     
    end
    
  
    
end
    

end


