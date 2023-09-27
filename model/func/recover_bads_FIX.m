function recover_bads_FIX(subIDtoFit, recovInfo, funInputFurther, sim_path)

path2save = strcat(sim_path, '/Fit', recovInfo.mod_str); 
mkdir(path2save)
des = recovInfo.exp_info;
acrR = recovInfo.exp_simcho; 
numSim = size(acrR, 2); 
options = bads('defaults');


if ( strcmp(recovInfo.mod_str , 'CFIX') )

    prior_mu = [0.5, 0     ];
    prior_var = [0.1,0.25   ];

    guessIn_LB   = [0.1,  -2 ];
    guessIn_UB   = [2.5,    2  ];

elseif  ( strcmp(recovInfo.mod_str , 'ZFIX') )
    prior_mu = [0.5  0  ];
    prior_var = [0.1  0 ];

    guessIn_LB   = [0.1 0];
    guessIn_UB   = [2.5  0];

    
end

LB   = guessIn_LB;
UB   = guessIn_UB;     

PLB = prior_mu -  sqrt(prior_var);
PUB = prior_mu + sqrt(prior_var);

NumRndIntials = 1; % can increase if you want multiple initial starting points
[xSet, ~] = lhsdesign_modified(NumRndIntials, PLB, PUB); 

options.UncertaintyHandling = true;
options.SpecifyTargetNoise = true;
options.MaxFunEvals = 10000;
options.CompletePoll = 'on';
options.NoiseFinalSamples = 1000; 


options_ibs = ibslike('defaults');
options_ibs.Nreps = 5;            %  % can increase this; Try and have a SD ~ 1 (and below 3) 
options_ibs.ReturnStd  = true;      % 2nd output is SD (not variance!)
                
options_ibs.NegLogLikeThreshold  = - length(find(acrR(:,1) ~= 0))*log(0.5);
options_ibs.Vectorized = false; 


                for iSim = 1: numSim
                    badsResults = [];x = []; fval=[] ; exitflag = []; output =[]; fInput = []; selR = zeros(1,NumRndIntials);
                    fInput = funInputFurther; 
                    fInput.designM = des(:,:, iSim); % Simulation virtual experiment
                    fInput.whichSim = iSim; 
                    fInput.mod_str = recovInfo.mod_str; %%%%%%%%%%% IMPORTANT
                    R = acrR(:, iSim); 
                    llfun = @(theta) ibslike(@fitting_fixedNN, theta, R, fInput.designM, options_ibs, fInput); 
                    
                    
                    for  irand = 1 : NumRndIntials
                        x0 = xSet(irand,:); % x0 =  [0.5, 0 ];
                        exitflag = 0; iter = 0; 
                        while exitflag == 0
                        [x,fval,exitflag,output] = bads(llfun,x0,LB,UB,PLB,PUB,[],options);
                                if exitflag ~= 1        % Retry fit if the previous one did not converge
                                [x,fval,exitflag,output] = bads(llfun,x0,LB,UB,PLB,PUB,[],options);

                                end

                            badsResults(irand).fitted_params{iter+1} = x;
                            badsResults(irand).fval(iter+1) = fval;
                            iter = iter + 1; 
                        end
                        selR(irand) = min(badsResults(irand).fval); 

                    end
                       
                    [~, min_ind] = min(selR);
                    storeResults = badsResults(min_ind); 
                    elapsed = toc; 
                    
                    parsave_best([strcat(path2save, '/Fit', recovInfo.mod_str, '_to','Sim', recovInfo.simulated,   num2str(subIDtoFit) )], subIDtoFit, storeResults)
                    fprintf( '[bads done with %s, Subj#%d]: took %d rounds of iterations',  recovInfo.mod_str, subIDtoFit, iter);
                    fprintf( ' : %.3f sec executed\n', elapsed );
                        
                      
                end
                
end
               