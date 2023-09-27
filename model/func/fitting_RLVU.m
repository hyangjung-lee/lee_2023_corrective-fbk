%% fitting script for model RLVU
function [R]= fitting_RLVU( theta, T , InputFurther)
whichSim = InputFurther.whichSim;  
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx; 
mod_str = InputFurther.mod_str; 

numCond =  InputFurther.numCond; 
repWindow = InputFurther.repWindow ;

 
numBlock_perRun=34;
numRun=InputFurther.numRun;
numTrial_perBlock=5;

tic;
 
learningRate =  theta(1); 
beliefNoiseSTD = theta(2);
z_mu_init = theta(3);
beta = theta(4); % inverse temperature parameter
init_V = theta(5);
extraRewardVal = InputFurther.extraRewardVal; 

lapse = InputFurther.lapseRate;

%%
 COND = unique(T(:,3)') ;

typeExtraReward = {'small','none', 'large'};			


    numTrial_perRun = numBlock_perRun*numTrial_perBlock;
    mR = []; 
 
    for iCond = 1:numCond
        
        iCond_typeER =  repmat({typeExtraReward{iCond}}, [numRun, numTrial_perRun]); 


        seqS = matT{subIDtoFit, iCond}.S(1:numRun, : );
        seqF = matT{subIDtoFit, iCond}.F(1:numRun, : ); % means the "CL" class variable;
        seqC = matT{subIDtoFit, iCond}.C(1:numRun, : , whichSim); 
        seqRT = matT{subIDtoFit, iCond}.RT(1:numRun,:); 
        seqC(seqRT > repWindow) = 0;
        lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 

        for iRun = 1 : size(seqS, 1)
           lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
           lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
        end
        wh_size =  size(seqS);  pC_vect = zeros(wh_size); 



        stimTrials = seqS;


        % set run numbers
        iterN = numRun;  % just one-shot
        trialN = length(stimTrials);

        % initialise variables, for speed
        pC = zeros(trialN,iterN);
        predictionError		= zeros(trialN,iterN);    
        V_L = NaN(trialN, iterN); V_R = NaN(trialN, iterN); 

     
         if (ismember( iCond, COND) )
                        
m_rnd = normrnd(0, beliefNoiseSTD, numRun, numTrial_perRun); 
            for iTrial =  1: max(T(:,2))
                extraReward = iCond_typeER{iRun, iTrial};

                % set contrast
                currentStim = stimTrials(:,iTrial)';
                % add sensory noise
                stim_withnoise  = currentStim + m_rnd(:, iTrial)';
                % calculate belief
                Belief_L = normcdf(z_mu_init,stim_withnoise,beliefNoiseSTD);
                Belief_R = 1 - Belief_L;
                if (iTrial == 1)
                    val_L = init_V*ones(size(currentStim)); 
                    val_R = init_V*ones(size(currentStim)) ; 

                else
                    vl = V_L(iTrial,:);  val_L = vl; % stored value of action L
                    vr = V_R(iTrial,:);  val_R = vr;  % stored value of action R
                end

ql = Belief_L .* val_L;
qr = Belief_R .* val_R;

%%%%%%%%%%%%%%%%%%%%%%%%% Decision Stage %%%%%%%%%%%%%%%%%%%%%%%%%       

% Softmax at action stage

maxQ = max([qr*beta; ql*beta]);
prob_stochL = exp( qr*beta - maxQ ) ./ ( exp( ql*beta - maxQ) + exp( qr*beta - maxQ) );
stochL = prob_stochL;
                
                
ch_tmp =  double( rand(1,numRun) < stochL );
ch_tmp(ch_tmp == 0)= -1;
pc = ch_tmp; 

    
pC(iTrial,:) = pc; 
pC(iTrial, find(seqC(:, iTrial)== 0) ) = 0; 
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
%  % calculate the prediction error, and update Q values
if strcmp(mod_str, 'RLVU')
    
    if (iTrial== 1)
       V_L(iTrial,:) = init_V;
       V_R(iTrial,:) = init_V;
       
    end
       V_L(iTrial + 1,:) =   V_L(iTrial,:);
       V_R(iTrial + 1,:) =    V_R(iTrial,:);
       
    
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




  
            

 
        
            end
    
         end

pC_vect = pC'; 
  
mR = [mR, pC_vect'];       

    
    end

vectR =    mR(:)  ;
R = vectR(T(:,1));

    

end


