function matModel = models_simulate(dataset,funInputFurther, numSim, theta_in, simulator_core)

iSub = funInputFurther.subindx;
imod = funInputFurther.imod;
mod_str = funInputFurther.mod_str; 

funInputFurther.data = dataset;
funInputFurther.numSim = numSim; 



    if (imod == 2 || imod == 3)
        funInputFurther.extraRewardVal = 0; 
    end
       
        tic;
        matModel = simulator_core(theta_in, [], funInputFurther);
        elapsed = toc; 

        fprintf( '[% s, subID#%d: sim done; %.3f sec]\n', mod_str, iSub, elapsed); 
        

    
end