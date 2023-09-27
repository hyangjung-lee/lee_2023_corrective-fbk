%% RLVU simulator
function [matModel]= simulator_RLVU( theta,~ , InputFurther)


 
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx; 
mod_str = InputFurther.mod_str; 
% imod = InputFurther.imod; 
% numCond =  InputFurther.numCond; 
numSim = InputFurther.numSim; 
 
numBlock_perRun   = 34;
numRun            = InputFurther.numRun;
numTrial_perBlock = 5;

    tic;
    
   

learningRate =  theta(1); % fixed at an initial value

beliefNoiseSTD = theta(2);


z_mu_init = theta(3);
 
 beta = theta(4); % inverse temperature parameter
 
extraRewardVal = InputFurther.extraRewardVal; 

 lapse = InputFurther.lapseRate;
%   lapse = 0; 
  init_V = theta(5); % 0.05

%% correction
  typeExtraReward = {'small','none', 'large'};			% o
%%
% COND = unique(T(:,3)') ;


    % Unbiased condition only
    numTrial_perRun=numBlock_perRun*numTrial_perBlock;
   sim_concat_m = []; 
for iSim = 1 : numSim    
    for iCond = 1:3
        
            iCond_typeER =  repmat({typeExtraReward{iCond}}, [10, 170]); 
            
%        for iRun = 1: numRun
 
            seqS = matT{subIDtoFit, iCond}.S;  vect_mu_c = zeros(size(seqS)); 
            seqF = matT{subIDtoFit, iCond}.F;
            seqC = matT{subIDtoFit, iCond}.C; 
          
lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 
            
            for iRun = 1 : size(seqS, 1)
               lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
               
               lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
            end
            pC_vect = zeros(size(seqS));  pL_vect_ev = zeros(size(seqS));
            confid = zeros(size(seqS)); 
% 
% seqS(seqS== 2) = 0.8;
% seqS(seqS== 1) = 0.4;
% seqS(seqS== 0) = 0;
% seqS(seqS== -1) = -0.4;
% seqS(seqS== -2) = -0.8;


            stimTrials = seqS;


            % set run numbers
            iterN = numRun;  % just one-shot
            iter = 1:iterN; 
            trialN = length(stimTrials);


            % initialise variables, for speed
            action						= cell(trialN,iterN);
            QL								= zeros(trialN,iterN);
            QR								= zeros(trialN,iterN);
            pC = zeros(trialN,iterN);
            predictionError		= zeros(trialN,iterN);
V_L = NaN(trialN,iterN); 
V_R = NaN(trialN,iterN);
     

     
%         if (ismember( iCond, COND) )
            
            % initalise Q values for each iteration
            QLL(1,iter) = 0;
            QRR(1,iter) = 0;
            QLR(1,iter) = 0;
            QRL(1,iter) = 0;
            
 m_rnd = normrnd(0, beliefNoiseSTD, numRun, numTrial_perRun)  ;         
            
            for iTrial = 1: numTrial_perRun
%                    iTrial = T(it);
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


                                            V_L(iTrial,:) = 1; % stored value of action L at 1st trial 
                                            V_R(iTrial,:) =  V_L(iTrial,:); % stored value of action R at 1st trial 

                                         
                                        else
                                            
                                            vl = V_L(iTrial,:);  val_L = vl; % stored value of action L
                                            vr = V_R(iTrial,:);  val_R = vr;  % stored value of action R
                                                                                      
                                        end
                                        
           
                                         
                                        %initialise Q values for this iteration
                                        ql = Belief_L.*val_L ; % expected value of action L
                                        qr = Belief_R.*val_R ;% expected value of action R
                                        
                                        QL(iTrial,iter) = Belief_L.*val_L + Belief_R .* QRL ; % expectedValue weighted by Belief_L and Belief_R ======> pL*V_L (in Lak, 2020)
                                        QR(iTrial,iter) = Belief_L.*QLR + Belief_R.*val_R; % pR*V_R in Lak (2020) 



% soft max normalization to prevent Inf or NaN

maxQ = max([qr*beta; ql*beta]); 
prob_stochL = exp( qr*beta - maxQ ) ./ ( exp( ql*beta - maxQ) + exp( qr*beta - maxQ) ); 
stochL = prob_stochL;


ch_tmp =  double( rand(1,numRun) < stochL ); 
ch_tmp(ch_tmp == 0)= -1;
pc = ch_tmp; 




% predictedvalue( pc == -1 ) = ql( pc == -1 );
% predictedvalue( pc == +1 ) = qr( pc == +1 ); 
% QC(iTrial, :) = predictedvalue;
pC(iTrial,:) = pc; 
                                        % chosen value Vc: value of chosen action
                                        vC(iTrial, pC(iTrial, :) == -1) = val_L(pC(iTrial, :) == -1);
                                        vC(iTrial, pC(iTrial, :) == +1) = val_R(pC(iTrial, :) == +1);
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
                                  
                                  Rd(iTrial, :) = reward ; 




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
                deltapredictionError(iTrial, iRun) = ( learningRate * pe );
            elseif (pC(iTrial, iRun) == +1)  % 'large'
                pe = reward(iRun) - V_R(iTrial+ 1, iRun) * Belief_R(iRun); 
                V_R(iTrial+1, iRun) =    V_R(iTrial, iRun) + learningRate * pe; 
                predictionError(iTrial, iRun) = pe;
                deltapredictionError(iTrial, iRun) =  (  learningRate * pe );
            end   








        end   
  
end
mat_reward(iTrial, :) = reward; 
mat_Belief_L(iTrial, :) = Belief_L; 
mat_mm(iTrial, :) = stim_withnoise;
  
            end

           
        

  pC_vect = pC'; 

            matModel{iCond}.C(:,:,iSim) = pC_vect;
            matModel{iCond}.F = seqF;
            matModel{iCond}.S = seqS;
            matModel{iCond}.RT = matT{subIDtoFit, iCond}.RT;
     
    end
    
  
    
end

end


