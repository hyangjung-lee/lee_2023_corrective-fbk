%% Boundary Likelihood P(m, CL | B)
function [Prob] = get_BoundLkld(z_range,  kalman_m,  I, sig_m, sig_mk, sig_S)


sig_M = sqrt( sig_mk^2 + sig_m^2 ); %% the sigma of p(m'|S) 
var_CDF = ( sig_M^2 * sig_S^2 ) / ( sig_M^2 + sig_S^2 );
sig_CDF = sqrt(var_CDF); 
           

    
term_pdf = normpdf(z_range, kalman_m, sqrt( sig_M^2 + sig_S^2 ));
term_cdf = NaN(size(term_pdf)); 
    
    
     
term_cdf(I == 1,:) =  1 - normcdf(z_range, (kalman_m(I == 1)*sig_S^2 + z_range*(sig_M)^2)/( sig_M^2 + sig_S^2 ), sig_CDF); 
term_cdf(I == -1,:) = normcdf(z_range, (kalman_m(I == -1)*sig_S^2 + z_range*(sig_M)^2)/( sig_M^2 + sig_S^2 ), sig_CDF); 
term_cdf(I == 0,:) = ones( length(find(I == 0)), length(z_range) ); 
       
          
            
Prob = term_cdf.*term_pdf;
    
         
   
end
