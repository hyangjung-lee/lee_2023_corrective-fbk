function [llkd, aicc, bic, prmbest] = load_gof_models(fitR, dtN)

pk      = size(fitR.fitted_params{1}, 2);
prmbest = fitR.fitted_params{1};

llkd    = -fitR.fval;  % LL
aicc    = -2*llkd + 2*pk + (2*pk*(pk+1))/(dtN-pk-1); % AICc
bic     = -2*llkd + pk*log(dtN) ; % BIC




end