function matModel = simulate_individual_lapgl(dataset,imod, mod_str, iSub, numSim, theta_in, lap_indv, simulator_core)
funInputFurther = []; x0=[];
subIDtoFit = iSub; 
funInputFurther.data = dataset;
funInputFurther.subindx = subIDtoFit;
funInputFurther.mod_str = mod_str; 
funInputFurther.imod = imod; 
funInputFurther.numCond = 3; 
funInputFurther.numRun = 10; 
funInputFurther.numSim = numSim; 
funInputFurther.lapseRate = lap_indv
 funInputFurther.extraRewardVal = 0; 

    if (imod == 2)
           funInputFurther.extraRewardVal = 0; 
        tic;
                matModel = simulator_core(theta_in, [], funInputFurther);
        elapsed = toc; 

            fprintf( '[ %s, subID#%d: sim done; %.3f sec]: learnRate=%.2f, , extraRV= %.2f,   beliefNoi=%.2f,  lapseRate =%.2f \n', ...
                    mod_str, iSub,elapsed,      theta_in(1),     funInputFurther.extraRewardVal ,       theta_in(2), funInputFurther.lapseRate);
                
                
    elseif (imod == 11)
        
        tic;
                matModel = simulator_core(theta_in, [], funInputFurther);
        elapsed = toc; 
            if (length(theta_in) == 1)
                   fprintf( '[ %s, subID#%d: sim done; %.3f sec]: mSig=%.2f,   lapseRate =%.2f \n', ...
                    mod_str, iSub,elapsed,      theta_in(1), funInputFurther.lapseRate) 
            else
            fprintf( '[ %s, subID#%d: sim done; %.3f sec]: mSig=%.2f, globalFixedBias=%.2f,  lapseRate =%.2f \n', ...
                    mod_str, iSub,elapsed,      theta_in(1),     theta_in(2),       theta_in(2), funInputFurther.lapseRate);
            end
    else

            tic;
                matModel = simulator_core(theta_in, [], funInputFurther);
        elapsed = toc; 

        
          fprintf( '[% s, subID#%d: sim done; %.3f sec]:mSig = %.2f,  z_mu_init= %.2f, lambda= %.2f,   sSig=%.2f, kSig= %.2f , noise_post = %.2f\n', ...
		mod_str, iSub, elapsed,                  theta_in(1),    theta_in(2),     theta_in(3),        theta_in(4),           theta_in(5),  theta_in(6));

    end
end