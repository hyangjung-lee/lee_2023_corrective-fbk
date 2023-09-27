%% RLVU simulator
function [matModel]= simulator_RLVU( theta,~ , InputFurther)


 
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx; 
mod_str = InputFurther.mod_str; 
% imod = InputFurther.imod; 
numCond =  InputFurther.numCond; 
numSim = InputFurther.numSim; 
 
numBlock_perRun   = 34;
numRun            = InputFurther.numRun;
numTrial_perBlock = 5;


learningRate =  theta(1); 
beliefNoiseSTD = theta(2);
z_mu_init = theta(3);
beta = theta(4); % inverse temperature parameter
init_V = theta(5); % 

extraRewardVal = InputFurther.extraRewardVal; 
lapse = InputFurther.lapseRate;

 

%% correction
typeExtraReward = {'small','none', 'large'};			



  
numTrial_perRun = numBlock_perRun*numTrial_perBlock;
 
for iSim = 1 : numSim    
    for iCond = 1:numCond
        
            iCond_typeER =  repmat({typeExtraReward{iCond}}, [numRun, numTrial_perRun]); 
 
            seqS = matT{subIDtoFit, iCond}.S(1:numRun, : );  
            seqF = matT{subIDtoFit, iCond}.F(1:numRun, : ); % means the "CL" class variable  
          
            lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS));         
            for iRun = 1 : size(seqS, 1)
               lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
               lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
            end



            stimTrials = seqS;


            % set run numbers
            iterN = numRun;  % just one-shot
            trialN = length(stimTrials);



            pC = zeros(trialN,iterN);
            predictionError		= zeros(trialN,iterN);
            V_L = NaN(trialN,iterN); V_R = NaN(trialN,iterN);
     

   
            
m_rnd = normrnd(0, beliefNoiseSTD, numRun, numTrial_perRun)  ;         
            
            for iTrial = 1: numTrial_perRun
                extraReward = iCond_typeER{iRun, iTrial};

                % set contrast
                currentStim = stimTrials(:,iTrial)';

                % add sensory noise
                stim_withnoise  = currentStim + m_rnd(:,iTrial)';


                % calculate belief
                    Belief_L = normcdf(z_mu_init,stim_withnoise,beliefNoiseSTD);
                    Belief_R = 1 - Belief_L;


                    if (iTrial == 1)

                        val_L = init_V*ones(size(currentStim)); 
                        val_R = init_V*ones(size(currentStim)); 


                    else

                        vl = V_L(iTrial,:);  val_L = vl; % stored value of action L
                        vr = V_R(iTrial,:);  val_R = vr;  % stored value of action R

                    end



%initialise Q values for this iteration
ql = Belief_L.*val_L ; % expected value of action L;  weighted by Belief_L and Belief_R ======> pL*V_L (in Lak, 2020)
qr = Belief_R.*val_R ;% expected value of action R;   weighted by Belief_L and Belief_R ======> pR*V_R  (in Lak, 2020)
               

% soft max normalization to prevent Inf or NaN
maxQ = max([qr*beta; ql*beta]); 
prob_stochL = exp( qr*beta - maxQ ) ./ ( exp( ql*beta - maxQ) + exp( qr*beta - maxQ) ); 
stochL = prob_stochL;


ch_tmp =  double( rand(1,numRun) < stochL ); 
ch_tmp(ch_tmp == 0)= -1;
pc = ch_tmp; 



pC(iTrial,:) = pc; 

%% whatever the QL & QR difference, lapse trials were made 

if (sum(lapse_idx(:,iTrial)) ~= 0)
    lap_trial_idx = lapse_idx(:,iTrial)==1; 
    pC(iTrial, lap_trial_idx) = lapse_val(lap_trial_idx,iTrial); 

    %%%% This part allows no update at all for the lapse trials
    Belief_L(lap_trial_idx) = 0;
    Belief_R(lap_trial_idx) = 0; 
end                      
                                  
reward = double( seqF(: , iTrial)' == pC(iTrial,:) ); 

if strcmp(extraReward, 'small')
    reward(seqF(: , iTrial)' == -1 &  pC(iTrial,:) == -1) =    reward(seqF(: , iTrial)' == -1 &  pC(iTrial,:) == -1) + extraRewardVal;
elseif strcmp(extraReward, 'large')
    reward(seqF(: , iTrial)' == +1 &  pC(iTrial,:) == +1) =    reward(seqF(: , iTrial)' == +1 &  pC(iTrial,:) == +1) + extraRewardVal;
end

if (strcmp(mod_str, 'RLVU'))
        if (iTrial == 1)
            V_L(iTrial,:) = init_V; 
            V_R(iTrial,:) = init_V; 
        end
            V_L(iTrial + 1, :) =   V_L(iTrial,:);
            V_R(iTrial + 1, :) =   V_R(iTrial,:);

        for iRun = 1 : numRun
            if (pC(iTrial, iRun) == -1)
                pe = reward(iRun) - V_L(iTrial+ 1, iRun) * Belief_L(iRun); 
                V_L(iTrial+1, iRun) =    V_L(iTrial, iRun) + learningRate * pe ; 
                predictionError(iTrial, iRun) = pe;
            elseif (pC(iTrial, iRun) == +1)  % 'large'
                pe = reward(iRun) - V_R(iTrial+ 1, iRun) * Belief_R(iRun); 
                V_R(iTrial+1, iRun) =    V_R(iTrial, iRun) + learningRate * pe; 
                predictionError(iTrial, iRun) = pe;
            end   



        end   
  
end
% mat_reward(iTrial, :) = reward; 
% mat_Belief_L(iTrial, :) = Belief_L; 
% mat_mm(iTrial, :) = stim_withnoise;
  
            end
        

  pC_vect = pC'; 

            matModel{iCond}.C(:,:,iSim) = pC_vect;
            matModel{iCond}.F = seqF;
            matModel{iCond}.S = seqS;
            matModel{iCond}.RT = matT{subIDtoFit, iCond}.RT;
     
    end
    
  
    
end

end


