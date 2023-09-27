%% Figure 1B
clear all; 

% color spec.
cmap_model(1,:) = [86 137 56]/255; cmap_model(2,:) = [160 63 116]/255;
% define axis
z_range = [-5:0.001:5]; x_axis = [0 2.5]; y_axis = [-0.2 0.2]; 


figure(11); clf; 
subplot(1,2,1); imodel = 2; 
    plot(x_axis,  [0 0 ], 'k--', 'linew', 2); hold on; plot( z_range, 0.3*( 1-normcdf( z_range, 0, 1) ), 'linew', 4, 'color', cmap_model(imodel,:)); box off
    xlim(x_axis) ; ylim(y_axis); set(gca,'ticklength', [0.03 0.03], 'xtick', [x_axis(1) x_axis(2)] , 'xticklabel', {'0','+\infty'}, 'linew', 2, 'ytick', [y_axis(1), 0, y_axis(2)], 'yticklabel',{'-','0', '+'}, 'fontsize', 20);
title('Value-updating')    

subplot(1,2,2); imodel = 1; 
    plot(x_axis,  [0 0 ], 'k--', 'linew', 2); hold on; plot( z_range, ( 1-normcdf( z_range, 0.3, 2) ) - 0.4, 'linew', 4, 'color', cmap_model(imodel,:)); box off
    xlim(x_axis) ; ylim(y_axis); set(gca,'ticklength', [0.03 0.03], 'xtick', [x_axis(1) x_axis(2)] , 'xticklabel', {'0','+\infty'}, 'linew', 2, 'ytick', [y_axis(1), 0, y_axis(2)], 'yticklabel',{'-','0', '+'}, 'fontsize', 20);
title('World-updating')    


%% PDF an example (conceptual)
figure(12); clf;
yobj = normpdf(z_range, 0, 1); 
plot(z_range,  yobj, 'k-','linew',2); hold on;  
area( z_range(z_range>0), yobj(z_range>0), 0, 'EdgeColor','none','FaceColor','k')
plot([0 0], [0 0.55],'k-','linew',4); axis off