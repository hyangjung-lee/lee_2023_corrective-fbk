clear all; close all; 
iscode = 1;

color_yellow = [0.9290 0.6940, 0.1250];
cmap_code{1} = color_yellow; 
cmap_code{2} = [255, 147, 0]/255; 
cmap_code{3} = [148, 33, 146]/255; 
cmap_code{4} = 0.8*[1 1 1]; 

z_range = [-5:.01:5]; 

%%%% Parameter set 
z_mu_init = 0;
z_sig_init = 1.5;
mSig = 0.5; 
sSig = 1.4; 
kSig = 1;
noise_post = 1;  
priSig = z_sig_init; 
%% example trial
iTrial = 1; 
cm = 0; 


m = 0.3071;
prevfbk = (+1);


mdashset = [0.1, 1.1, 2.1]; 
%% 

E_c_dS = cm; E_x_dS = m; 

choice =  double(  E_x_dS - E_c_dS > 0 );
choice(choice == 0) = -1 ;

prevdec = choice;


if (prevfbk == -1)
    if (choice == -1)
        Fgiven = 'correct';
    elseif (choice == +1)
        Fgiven = 'incorrect';
    end
elseif (prevfbk == +1)
    if (choice == -1)
        Fgiven = 'incorrect';
    elseif (choice == +1)
        Fgiven = 'correct';
    end
end

for iK = 1 : length(mdashset)


    if prevdec > 0
        kalman_m = mdashset(iK);  

    else
        kalman_m = - mdashset(iK);  

    end



    fd_match = prevfbk == prevdec; 
    I(fd_match) = prevfbk(fd_match); I(~fd_match) =  prevfbk(~fd_match);


    sig_kmL = sqrt( kSig^2 + mSig^2 ); %% the sigma of p(m'|S) 
    var_CDF = (sSig^2 * sig_kmL^2) / (sSig^2 + sig_kmL^2);
    sig_CDF = sqrt(var_CDF); 

    term_pdf = normpdf(z_range, kalman_m, sqrt( sSig^2 + kSig^2 + mSig^2 ));

    term_cdf = NaN(size(term_pdf)); 

    if (I == -1)
        term_cdf = normcdf(z_range, (kalman_m*sSig^2 + z_range*(sig_kmL)^2)/(sSig^2 + sig_kmL^2), sig_CDF);
    elseif (I == +1 )
        term_cdf = 1 - normcdf(z_range, (kalman_m*sSig^2 + z_range*(sig_kmL)^2)/(sSig^2 + sig_kmL^2), sig_CDF); 
    else
        term_cdf = ones( 1, length(z_range) ); 
    end


    Prob = term_cdf.*term_pdf;
    sI = sum(Prob,2);  norm_Prob = Prob./sI; 
    like_mu = sum(norm_Prob .* z_range, 2);%%%%%%%% the mean of distribution ( mu_Prob : Normalized posterior  ) 

    pri =  normpdf(z_range, E_c_dS, priSig); 



    postProb = pri .* Prob ; % step 2: mu posterior
    sI = sum(postProb,2);  norm_postProb = postProb./sI; 
    cfw_mu = sum(norm_postProb .* z_range, 2);%%%%%%%% the mean of distribution ( mu_Prob : Normalized posterior  ) 
    is_post = sqrt( sum(( (z_range -  cfw_mu ).^2).*norm_postProb, 2) ); %%%%%%%% the s.d of distribution (mu_Prob: Normalized posterior)


    % post-process (2) adding inclination to long-term prior 

    lambda = ( z_sig_init^2)./(z_sig_init^2 + is_post.^2);
    lambda( z_sig_init == 0 ) = 0; % when prior is really sharp 
    lambda( z_sig_init == Inf ) = 1; % when prior is really broad 



    vect_mu_c(:,iTrial) = cfw_mu .* lambda + z_mu_init .* (1-lambda); %% imu_pre = la*imu_post + (1-la) * z_mu_init



    figure(iK+110);   set(gcf,'position', [73 136 274 626]);

    if (abs(kalman_m) == 0.1)
        iscode = 1;  
    elseif (abs(kalman_m) == 1.1)
        iscode = 2; 
    elseif (abs(kalman_m) == 2.1)
        iscode = 3;   
    else
        iscode = 4; 
    end

    subplot(3,1,1);
        yobj = term_pdf./sum(term_pdf); yobj_mu = sum(yobj .* z_range, 2);
        plot(z_range, yobj , 'color', cmap_code{iscode}); hold on;  
        yd = yobj; area( z_range, yd, 0, 'EdgeColor','none','FaceColor',cmap_code{iscode}, 'facealpha',0.6);

        plot([0 0 ], [0 0.5], 'k--','linew',1);
        plot([yobj_mu yobj_mu ], [0 0.5], 'k-','linew',2); ax1 = gca; ax1.YAxis.Visible = 'off'; ylim([0 0.0032]);xlim([-3 3]); 
    box off; set(gca,'fontsize', 18, 'xtick', [-2:2], 'ticklength', [0.02 0.03]); 

    subplot(3,1,2);

        plot(z_range,term_cdf, 'color', cmap_code{iscode}, 'linew', 2); hold on;  
        plot([0 0 ], [0 1], 'k--','linew',1); ax1 = gca; ax1.YAxis.Visible = 'off';xlim([-3 3]);
    box off; set(gca,'fontsize', 18, 'xtick', [-2:2], 'ticklength', [0.02 0.03]);


    subplot(3,1,3);
        yobj = Prob./sum(Prob); yobj_mu = sum(yobj .* z_range, 2);
        plot(z_range, yobj , 'color', cmap_code{iscode}); hold on;  

        yd = yobj; area( z_range, yd, 0, 'EdgeColor','none','FaceColor',cmap_code{iscode}, 'facealpha',0.6);

        plot([0 0 ], [0 0.5], 'k--','linew',1);
        plot([yobj_mu yobj_mu ], [0 0.5], 'k-','linew',2); ax1 = gca; ax1.YAxis.Visible = 'off'; ylim([0 0.0032]);xlim([-3 3])
    box off; set(gca,'fontsize', 18, 'xtick', [-2:2], 'ticklength', [0.02 0.03]); 

end