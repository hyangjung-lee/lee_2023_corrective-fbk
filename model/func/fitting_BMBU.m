%% fitting script for model BMBU 
function [R]= fitting_BMBU( theta, T , InputFurther)

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

mSig =  theta(1); 
z_mu_init = theta(2);
z_sig_init = theta(3); 
sSig = theta(4);
kSig =  theta(5); %% noise m -> mm (memory noise always bigger than mSig) 
noise_post =  theta(6); % diffusion param

lapse = InputFurther.lapseRate;

%%
COND = unique(T(:,3)') ;

mu0 = 0; ww = 10;
ran_mu = 10;resol_mu = 200/ww; % resolution of mu
z_range = linspace(mu0 - ran_mu, mu0 + ran_mu, resol_mu); len_ax = length(z_range); 
Prob = nan(1,len_ax);



    numTrial_perRun = numBlock_perRun*numTrial_perBlock;

    mR = [];
   
    for iCond = 1:numCond
        seqS = matT{subIDtoFit, iCond}.S(1:numRun, : );
        lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 
        for iRun = 1 : size(seqS, 1)
           lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse;   
           lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
        end
        seqF = matT{subIDtoFit, iCond}.F(1:numRun, : ); % means the "CL" class variable
        seqC = matT{subIDtoFit, iCond}.C(1:numRun, : , whichSim);
        seqRT = matT{subIDtoFit, iCond}.RT(1:numRun,:); 
        seqC(seqRT > repWindow) = 0; % choice not count (miss)
        I = NaN(numRun, 1); mu_Prob = NaN(numRun, length(z_range)); postProb =  NaN(numRun, length(z_range));
        wh_size =  size(seqS); vect_mu_c = zeros(wh_size); pC_vect = zeros(wh_size); vect_SD = zeros(wh_size);
         
m_rnd = normrnd(0, mSig, numRun, numTrial_perRun);        
km_rnd = normrnd(0, kSig, numRun, numTrial_perRun);        

     
        if (ismember( iCond, COND) )
            for iTrial = 1: max(T(:,2))
                if(iTrial == 1)
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
                    priUmpSig(priUmpSig == 0) = 0.1; % jump out 
cm   =   preC;
x_mu = seqS(:,iTrial);
m    = x_mu + m_rnd(:, iTrial);
kalman_m = m + km_rnd(:,iTrial);

                  
%%%%%%%%%%%%%%%%%%%%%%%%% Decision Stage %%%%%%%%%%%%%%%%%%%% algthm. equiv.            
pC_vect( m - cm > 0 , iTrial) = 1;
pC_vect( m - cm < 0  , iTrial) = -1;
pC_vect( m == cm  , iTrial) = 0 ;
                   
%%%% let model do the same thing for overdue response time
pC_vect( seqC(:, iTrial) == 0 , iTrial) = 0;


if (sum(lapse_idx(:,iTrial)) ~= 0)
    lap_trial_idx = lapse_idx(:,iTrial)==1; 
    pC_vect(lap_trial_idx, iTrial) = lapse_val(lap_trial_idx,iTrial); 
end

     prevfbk = seqF(:,iTrial);
     prevdec = pC_vect(:,iTrial);%%%% this is modified from seqC to pC_vect (after Slide 20) 

     

    fd_match = prevfbk == prevdec; 
    I(fd_match) = prevfbk(fd_match); I(~fd_match) =  prevfbk(~fd_match);
    updSig = priUmpSig;
 


%%%  Boundary Likelihood
 
    Prob = get_BoundLkld( z_range, kalman_m, I, mSig, kSig, sSig);

%%% Boundary Prior      
              
    pri =  normpdf(z_range, cm, updSig);
              


               
                if (sum(lapse_idx(:,iTrial)) ~= 0)
                    postProb(lap_trial_idx,:) = pri(lap_trial_idx,:) ; 
                    postProb(~lap_trial_idx,:) = pri(~lap_trial_idx,:) .* Prob(~lap_trial_idx,:); % step 2: mu posterior

                else

                    postProb = pri .* Prob ; % step 2: mu posterior

                end

                    sI = sum(postProb,2);  norm_postProb = postProb./sI; 
                    val_idx = ( sI ~= len_ax ) ; non_val = find( sI == len_ax); 
                    mu_Prob(val_idx,:) = norm_postProb(val_idx,:) ;

                     if (~isempty(non_val))
                        mu_Prob(non_val,:) = normpdf(z_range, z_mu_init, kSig)./sum(normpdf(z_range, z_mu_init, kSig)); 
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
 



            end

           
        end



mR = [mR, pC_vect']; 




    end
    
    
vectR =    mR(:)  ;



R = vectR(T(:,1));

end


