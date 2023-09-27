function [aftFit, stand_aftFit] = fit_psych_PSE(data,  options,  iboot ) 


nTri = data(:,3);
yyA = data(:,2)./data(:,3);
nanresult = NaN(5,1);
if (size(data, 1) < 3)
     aftFit = nanresult; 
     stand_aftFit = nanresult; 
else
     if ( isempty(iboot) || iboot == 1 )
         result = psignifit(data,options);
         aftFit = result.Fit;
     else
 
        nCho_bt = binornd( nTri, yyA );
        
        data(:,2) = nCho_bt;
        result = psignifit(data,options);
        aftFit = result.Fit;
     end
       stand_aftFit = getStandardParameters(result.Fit, 'norm', 0.05); 

       
%      clf
     result.data = data;
%    plotPsych(result); drawnow;

end


end
