%% Plotting # of trials
clear all; 
cd ../
path_str.currpath   = pwd;

which_toi =  +1 

% define 1-D axis of a matrix to plot (toi)
axis_id = [1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20]; 

which_data = 1
    if (which_data == 1)
        path_str.data       = [path_str.currpath, '/data']; 
        load([path_str.data, '/Episode_freq_OBS.mat'])
    end
    if (which_toi == -1)
        CUB = toiminusone_Pcub;
    elseif (which_toi == +1)
        CUB = toiplusone_Pcub;
    end



for iSub = 1 : 30
    indChunk = []; 
    for id = 1 : length(axis_id)

        if (isnan(axis_id(id)))
            indChunk(:, id) = NaN(size(axis_id));

        else
            idch = axis_id(id); 

                if (idch == 1) 
                    [norm_cube] = ( CUB.extPure{2}.SMALL(:,:,:, iSub) );% XS -> small -> cor(veridical) (when this toi clamped)
                elseif (idch == 2)
                    [norm_cube] = ( CUB.norPure{2}.SMALL(:,:,:, iSub) );% S -> small -> cor(veridical)
                elseif (idch == 3)
                    [norm_cube] = ( CUB.extPure{6}.SMALL(:,:,:, iSub) );% M -> small -> cor
                elseif (idch == 4) 
                    [norm_cube] = ( CUB.norPure{4}.SMALL(:,:,:, iSub) );% L -> small -> cor
                elseif (idch == 5)
                    [norm_cube] = ( CUB.extPure{4}.SMALL(:,:,:, iSub) );% XL -> small -> cor


                elseif (idch == 6)
                    [norm_cube] = ( CUB.extPure{4}.LARGE(:,:,:, iSub) );% XS -> large -> cor
                elseif (idch == 7)
                    [norm_cube] = ( CUB.norPure{4}.LARGE(:,:,:, iSub) ); % S -> large -> cor
                elseif (idch == 8)
                    [norm_cube] = ( CUB.extPure{6}.LARGE(:,:,:, iSub) ); % M -> large -> cor
                elseif (idch == 9)   
                    [norm_cube] = ( CUB.norPure{2}.LARGE(:,:,:, iSub) ); % L -> large -> cor (veridical)
                elseif (idch == 10)
                    [norm_cube] = ( CUB.extPure{2}.LARGE(:,:,:, iSub) );% XL -> large -> cor(veridical)

                elseif (idch == 11)
                    [norm_cube] = ( CUB.extPure{3}.SMALL(:,:,:, iSub) ); % XS -> small -> incor
                elseif (idch == 12)
                    [norm_cube] = ( CUB.norPure{3}.SMALL(:,:,:, iSub) );% S -> small -> incor
                elseif (idch == 13)
                    [norm_cube] = ( CUB.extPure{7}.SMALL(:,:,:, iSub) );% M -> small -> incor
                elseif (idch == 14) 
                    [norm_cube] = ( CUB.norPure{5}.SMALL(:,:,:, iSub) ); % L -> small -> incor
                elseif (idch == 15)
                    [norm_cube] = ( CUB.extPure{5}.SMALL(:,:,:, iSub) );% XL -> small -> incor


                elseif (idch == 16)
                    [norm_cube] = ( CUB.extPure{5}.LARGE(:,:,:, iSub) );% XS -> large -> incor
                elseif (idch == 17)
                    [norm_cube] = ( CUB.norPure{5}.LARGE(:,:,:, iSub) );% S -> large -> incor
                elseif (idch == 18)
                    [norm_cube] = ( CUB.extPure{7}.LARGE(:,:,:, iSub) );% M -> large -> incor
                elseif (idch == 19)   
                    [norm_cube] = ( CUB.norPure{3}.LARGE(:,:,:, iSub) );% L -> large -> incor
                elseif (idch == 20)
                    [norm_cube] = ( CUB.extPure{3}.LARGE(:,:,:, iSub) );% XL -> large -> incor



                end
                
                
                clear episodefreq
                for iS = 1:5 % XS to XL
                iCh = 1; % small
                iFb = 2; % correct
                episodefreq(iS) = norm_cube(iS , iCh, iFb);
                end
                chunk1 = episodefreq; 
                short_chunk1 = [ sum( chunk1(1:2) ),  chunk1(3), sum(chunk1(4:5))]; 
                for iS = 1:5 % XS to XL
                iCh = 2; % large
                iFb = 2; % correct
                episodefreq(iS) = norm_cube(iS , iCh, iFb);
                end
                chunk2 = episodefreq; 


                for iS = 1:5 % XS to XL
                iCh = 1; % small
                iFb = 1; % incorrect
                episodefreq(iS) = norm_cube(iS , iCh, iFb);
                end
                chunk3 = episodefreq; 


                for iS = 1:5 % XS to XL
                iCh = 2; % large
                iFb = 1; % incorrect
                episodefreq(iS) = norm_cube(iS , iCh, iFb);
                end
                chunk4 = episodefreq; 


            stack_idx = find( axis_id == axis_id(id) ); 
            indChunk(:,stack_idx ) = [chunk1,  chunk2,  chunk3, chunk4];
        end

    end

    if (which_data == 1)
%         title('Humans')
         ALL_indvOBS(:,:,:,iSub) = indChunk; 
   
    end


 
end





raw_h = squeeze( sum(ALL_indvOBS, 1) )'; 
axis_id = [1,2,3,4,5, NaN,6,7,8,9,10, NaN, 11,12,13,14,15, NaN, 16,17,18,19,20]; 
ALL_raw_h = []; 
for id = 1 : size(axis_id, 2)
    if (isnan(axis_id(id)))
        ALL_raw_h(:, id) = NaN(30,1); 
        
    else
        j = axis_id(id); 
    
        
        ALL_raw_h(:, id) =       raw_h(:,j) ; 
        
    end
end


% %% Fig 7D
code_color = [0 0 0 ]; 
imodel = 3; 
figure(11+11);clf

cond = [1, 11];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [  1 1 1]; h.EdgeColor = code_color;  h.LineWidth = 1.5;  hold on; 
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '-','color', code_color , 'linew', 1.5); 
end

cond = [2,3, 9,10];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [  1 1 1]; h.EdgeColor =  code_color; h.LineWidth = 1.5;
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '--','color', code_color, 'linew', 1.5); 

end
cond = [4,5, 7,8];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [1 1 1]; h.EdgeColor =  code_color;h.LineWidth = 2; h.LineStyle = ':';
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '-','color', code_color, 'linew', 1); 
end


cond = [17, 19];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [  1 1 1]; h.EdgeColor = code_color;  h.LineWidth = 1.5;  hold on; 
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '-','color', code_color , 'linew', 1.5); 
end

cond = [15,16, 20,21];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [  1 1 1]; h.EdgeColor =  code_color; h.LineWidth = 1.5; 
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '--','color', code_color, 'linew', 1.5); 

end
cond = [13,14, 22, 23];
for j = 1: length(cond) 
    h=bar(cond(j), mean(ALL_raw_h(:,cond(j)))); h.FaceColor = [1 1 1]; h.EdgeColor = code_color; h.LineWidth = 2; h.LineStyle = ':';
    errorbar(cond(j),  mean(ALL_raw_h(:,cond(j))), std(ALL_raw_h(:,cond(j) ) ), '-','color', code_color, 'linew', 1); 
end


set(gca,'xtick', [], 'fontsize', 20); 
ylabel('Number of trials'); box off;  