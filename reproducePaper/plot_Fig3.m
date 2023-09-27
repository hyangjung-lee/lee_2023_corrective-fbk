%% Figure 3
% Figure 3A,B,C,D,E

clear all; close all; 
%% PDF measurement example (for Figure 3A)
z_range = [-5:0.001:5]; 

yd = normpdf( z_range, 0, 1); 
crit_val = -0.3;

figure(31); clf;  
plot(z_range, yd, 'color', 'k', 'linew',2 ); hold on;  
   
    ind = (yd> 0)& z_range< crit_val;
    yd(~ind) = NaN; 
    % coloring the left-side of the PDF 
    area( z_range, yd, 0, 'EdgeColor','none','FaceColor','none', 'facealpha',0.8);

    % coloring the right-side of the PDF
    yd = normpdf( z_range, 0, 1); 
    ind = (yd> 0)& z_range > crit_val;
   
    yd(~ind) = NaN; 
    area( z_range, yd, 0, 'EdgeColor','none','FaceColor','k', 'facealpha',0.8);

plot([crit_val crit_val], [0 0.4], 'k-','linew',2); ylabel('Probability');


%% for Figure 3 B,C Psychometric curves (conceptual) 
z_range = linspace(-10, 10, 2000);
color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];

figure(32);clf; 

plot(z_range, normcdf( z_range, -1, 1), 'color', color_yellow, 'linew', 3); xlim([-5 5]); hold on; 
plot([-5 5], [0.5 0.5], 'k--'); plot([0 0], [0 1], 'k--');

plot(z_range, normcdf( z_range, 0, 1), 'color', 'k', 'linew', 3); 
plot(z_range, normcdf( z_range, +1, 1), 'color', color_blue, 'linew', 3);

%% Figure 3D (conceptual) 
figure(33);clf; 
clear xx yy ffx ffy;
ffx = z_range( find ( z_range < 0) );
xx = [ffx; ffx];
yy = [ffx; 0*ffx];

s1 = surf(xx,yy, xx*0,'alphadata', fliplr(xx), 'facecolor', color_yellow, 'edgecolor','none', 'facealpha', 'interp'); hold on;
clear xx yy ffx ffy;

ffx = z_range( find ( z_range > 0) );
xx = [ffx; ffx];
yy = [ffx; 0*ffx];

s1 = surf(xx,yy, xx*0,'alphadata', xx, 'facecolor', color_blue, 'edgecolor','none', 'facealpha', 'interp' );
ylim([-10 10]); xlim([-10 10]); view(2)

%% Figure 3E (conceptual) 

figure(34);clf; 
clear xx yy ffx ffy;
ffx = z_range( find ( z_range < 10) );
xx = [ffx; ffx];
yy = [0.5*ffx-5; 0*ffx];

s1 = surf(xx,yy, xx*0,'alphadata', fliplr(xx), 'facecolor', color_yellow, 'edgecolor','none', 'facealpha', 'interp');
ylim([-10 10]); xlim([-10 10]); view(2)