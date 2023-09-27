%% BMBU simulator 
function [matModel]= simulator_BMBU( theta, ~ , InputFurther)
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx;  
mod_str = InputFurther.mod_str; 
% imod = InputFurther.imod; 

numCond =  InputFurther.numCond; 
numSim = InputFurther.numSim; 
 
numBlock_perRun   = 34;
numRun            = InputFurther.numRun;
numTrial_perBlock = 5;

 
mSig =  theta(1); 
 
z_mu_init = theta(2);
z_sig_init = theta(3); 

     
sSig= theta(4);

kSig = theta(5); %% noise m -> mm
noise_post =  theta(6); % diffusion param
 
 
lapse = InputFurther.lapseRate;

mu0 = 0; ww = 10;
ran_mu = 10;resol_mu = 200/ww; % resolution of mu
z_range = linspace(mu0 - ran_mu, mu0 + ran_mu, resol_mu); len_ax = length(z_range); 
Prob = nan(1,len_ax);


numTrial_perRun = numBlock_perRun*numTrial_perBlock;
 
for iSim = 1 : numSim    
   
    for iCond = 1: numCond


        seqS = matT{subIDtoFit, iCond}.S(1:numRun, : );

        lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 
        for iRun = 1 : size(seqS, 1)
           lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
           lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
        end

            
        seqF = matT{subIDtoFit, iCond}.F(1:numRun, : ); % means the "CL" class variable
        I = NaN(numRun, 1); mu_Prob = NaN(numRun, length(z_range)); 
          
        wh_size =  size(seqS); vect_mu_c = zeros(wh_size);pC_vect = zeros(wh_size);  vect_SD = zeros(wh_size);
%       vect_kalman_m =  zeros(wh_size);vect_crit  =  zeros(wh_size);vect_decm =  zeros(wh_size);


     
m_rnd = normrnd(0, mSig, numRun, numTrial_perRun);        
km_rnd = normrnd(0, kSig, numRun, numTrial_perRun);        

            for iTrial = 1:numTrial_perRun
                if(iTrial==1)
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        preC = repmat(z_mu_init, numRun, 1) ;  
                        priUmpSig = repmat( z_sig_init, numRun, 1) ;  
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                % updating fc based on stimulus, choice & feedback

                        preC = vect_mu_c(:,iTrial-1);
                        priUmpSig = vect_SD(:,iTrial - 1)  ;  
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end
                        priUmpSig(priUmpSig == 0) = 0.1; 


                cm              =   preC;


x_mu = seqS(:,iTrial);
m = x_mu + m_rnd(:, iTrial);
kalman_m = m + km_rnd(:,iTrial);




                  
%%%%%%%%%%%%%%%%%%%%%%%%% Decision Stage %%%%%%%%%%%%%%%%%%%%%%%%% algthm. equiv.       
pC_vect( m - cm > 0 , iTrial) = 1;
pC_vect( m - cm < 0  , iTrial) = -1;
pC_vect( m == cm  , iTrial) = 0 ;

            
 



if (sum(lapse_idx(:,iTrial)) ~= 0)
    lap_trial_idx = lapse_idx(:,iTrial)==1; 
    pC_vect(lap_trial_idx, iTrial) = lapse_val(lap_trial_idx,iTrial); 
    
end
 
     prevfbk = seqF(:,iTrial);
     prevdec = pC_vect(:,iTrial);%%%% this is modified from seqC to pC_vect (after Slide 20) 

if   (strcmp(mod_str, 'BMBU'))
   
    fd_match = prevfbk == prevdec; 
    I(fd_match) = prevfbk(fd_match); I(~fd_match) =  prevfbk(~fd_match);
    updSig = priUmpSig;
    Prob = get_BoundLkld( z_range, kalman_m, I, mSig, kSig, sSig);
end
 
 
 
     
 
%%% Boundary Prior     
    pri =  normpdf(z_range, cm, updSig); 
   
   
  


        if (sum(lapse_idx(:,iTrial)) ~= 0) %%% if any case, lapse trials are contained 
          
            postProb(lap_trial_idx,:) = pri(lap_trial_idx,:) ; 
            postProb(~lap_trial_idx,:) = pri(~lap_trial_idx,:) .* Prob(~lap_trial_idx,:); % step 2: mu posterior

        else % if there's no lapse in the given trial 

          postProb = pri .* Prob ; % step 2: mu posterior

        end
        sI = sum(postProb,2);  norm_postProb = postProb./sI; 

        val_idx = ( sI ~= len_ax ) ; 
        mu_Prob(val_idx,:) = norm_postProb(val_idx,:) ;
        non_val = find( sI == len_ax ); % if prior x likelihood = [1 ..... 1] -> go back to initial belief  
        if (~isempty(non_val))

            mu_Prob(non_val,:) = normpdf(z_range, z_mu_init, z_sig_init)./sum(normpdf(z_range, z_mu_init, z_sig_init)); 
        end
         
         


        cfw_mu = sum(mu_Prob .* z_range, 2);%%%%%%%% the mean of distribution ( mu_Prob : Normalized posterior  ) 
        is_post = sqrt( sum(( (z_range -  cfw_mu ).^2).*mu_Prob, 2) ); %%%%%%%% the s.d of distribution (mu_Prob: Normalized posterior)
 

% post-process:  inclination to long-term prior 
lambda = ( z_sig_init^2)./(z_sig_init^2 + is_post.^2);
lambda( z_sig_init == 0 ) = 0; % when prior is really sharp 
lambda( z_sig_init == Inf ) = 1; % when prior is really broad 



vect_mu_c(:,iTrial) = cfw_mu .* lambda + z_mu_init .* (1-lambda); %% imu_pre = la*imu_post + (1-la) * z_mu_init

diffus_term = noise_post; 
 
 
 
% post-process: + diffusion process
vect_SD(:,iTrial) = sqrt( (lambda .* (is_post.^2)) + diffus_term.^2 ) ;
 


% vect_kalman_m(:,iTrial) = kalman_m; 
% vect_crit(:, iTrial) = preC; 
% vect_decm(:, iTrial) = m; 
            end

           
        

%%%% variables out

        matModel{iCond}.C(:,:,iSim) = pC_vect;
        matModel{iCond}.S = seqS;

        matModel{iCond}.RT = matT{subIDtoFit, iCond}.RT;
        matModel{iCond}.F = seqF;

    end
end
    



end


