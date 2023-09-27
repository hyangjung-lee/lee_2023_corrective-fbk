%% fitting script for Fixed (CFIX ZFIX)
function [R]= fitting_fixedNN( theta, T , InputFurther)


whichSim = InputFurther.whichSim; 
matT = InputFurther.data;
subIDtoFit = InputFurther.subindx; 
mod_str = InputFurther.mod_str; 

numCond =  InputFurther.numCond; 
repWindow = InputFurther.repWindow ;
 
numBlock_perRun   = 34;
numRun            = InputFurther.numRun;
numTrial_perBlock = 5;

tic;  

mSig =  theta(1); 
z_mu_init = theta(2);


lapse = InputFurther.lapseRate;

%%
COND = unique(T(:,3)') ;

numTrial_perRun=numBlock_perRun*numTrial_perBlock;

mR = [];
   
    for iCond = 1:numCond
            seqS = matT{subIDtoFit, iCond}.S(1:numRun, : ); wh_size =  size(seqS); pC_vect = zeros(wh_size); 
            lapse_val = nan(size(seqS));   lapse_idx = nan(size(seqS)); 
            for iRun = 1 : size(seqS, 1)
               lapse_idx(iRun,:) = rand(size(seqS(iRun,:))) < lapse; 
               lapse_val(iRun,lapse_idx(iRun,:)== 1) = randi(2,[sum(lapse_idx(iRun,:)),1])*2-3; 
            end
            % seqF = matT{subIDtoFit, iCond}.F(1:numRun, : );
            seqC = matT{subIDtoFit, iCond}.C(1:numRun, : , whichSim);
            seqRT = matT{subIDtoFit, iCond}.RT(1:numRun,:); 
            seqC(seqRT > repWindow) = 0;
      
m_rnd = normrnd(0, mSig, numRun, numTrial_perRun);        
   

     
        if (ismember( iCond, COND) )
            for iTrial = 1: max(T(:,2))


               


x_mu = seqS(:,iTrial);
m = x_mu + m_rnd(:, iTrial);

 

cm = z_mu_init*ones(size(x_mu));



                  
 %%%%%%%%%%%%%%%%%%%%%%%%% Decision Stage %%%%%%%%%%%%%%%%%%%%%%%%%       
pC_vect( m - cm > 0 , iTrial) = 1;
pC_vect( m - cm < 0  , iTrial) = -1;
pC_vect( m == cm  , iTrial) = 0 ;

                   
%%%% let model do the same thing for overdue response time
pC_vect( seqC(:, iTrial) == 0 , iTrial) = 0;


if (sum(lapse_idx(:,iTrial)) ~= 0)
    lap_trial_idx = lapse_idx(:,iTrial)==1; 
    pC_vect(lap_trial_idx, iTrial) = lapse_val(lap_trial_idx,iTrial); 
end

            end

           
        end



mR = [mR, pC_vect']; 




    end
    
    
vectR =    mR(:)  ;



R = vectR(T(:,1));

end


