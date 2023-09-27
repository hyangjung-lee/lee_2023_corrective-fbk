clear all; close all;
cd ../

path_str.currpath   = pwd;
addpath([path_str.currpath, '/analysis_func'])
color_blue = [0.3010, 0.7450, 0.9330];
color_yellow = [0.9290 0.6940, 0.1250];
xxd = [-2:2]; Xxd = [ones(size(xxd')), xxd']; linealph = 1; 


%% plotting Fig S6
%% %% Get Population-level Episode Frequency Maps
which_toi = -1; % case 'Retrospective'
which_toi = +1; % case 'Prospective'


% define 1-D axis of a matrix to plot
axis_id = [1,2,3,4,5, NaN,6,7,8,9,10, NaN, 11,12,13,14,15, NaN, 16,17,18,19,20]; 
nStim   =  5; 
for which_data = 1:3
    if (which_data == 1)
        path_str.data       = [path_str.currpath, '/data']; 
        load([path_str.data, '/Episode_freq_OBS.mat'])
    elseif (which_data == 2)
        path_str.data       = [path_str.currpath, '/results_expost']; 
        load([path_str.data,'/BMBU/Episode_freq_Wmodel_verS445_FF.mat'])
    elseif (which_data == 3)
        path_str.data       = [path_str.currpath, '/results_expost']; 
        load([path_str.data, '/RLVU/Episode_freq_Vmodel_verS_softbias_plus_SSpercpt.mat'])
    end
    if (which_toi == -1)
        CUB = toiminusone_Pcub;
    elseif (which_toi == +1)
        CUB = toiplusone_Pcub;
    end


    indChunk = NaN(size(axis_id, 2), size(axis_id, 2)); % 23-by-23 matrix for colored maps in S6 Fig 
                    

    for id = 1 : length(axis_id)
        if (isnan(axis_id(id)))
            indChunk(:, id) = NaN(size(axis_id, 2), 1);
        else
            idch = axis_id(id); 
                % Episode type I: C = small and F = correct
                if (idch == 1)
                    [norm_cube] = get_normCube( CUB.extPure{2}.SMALL );% XS -> small -> cor(veridical) (when this toi clamped)
                elseif (idch == 2)
                    [norm_cube] = get_normCube( CUB.norPure{2}.SMALL );% S -> small -> cor(veridical)
                elseif (idch == 3)
                    [norm_cube] = get_normCube( CUB.extPure{6}.SMALL );% M -> small -> cor
                elseif (idch == 4) 
                    [norm_cube] = get_normCube( CUB.norPure{4}.SMALL );% L -> small -> cor
                elseif (idch == 5)
                    [norm_cube] = get_normCube( CUB.extPure{4}.SMALL );% XL -> small -> cor

                % Episode type II: C = large and F = correct
                elseif (idch == 6)
                    [norm_cube] = get_normCube( CUB.extPure{4}.LARGE );% XS -> large -> cor
                elseif (idch == 7)
                    [norm_cube] = get_normCube( CUB.norPure{4}.LARGE ); % S -> large -> cor
                elseif (idch == 8)
                    [norm_cube] = get_normCube( CUB.extPure{6}.LARGE ); % M -> large -> cor
                elseif (idch == 9)   
                    [norm_cube] = get_normCube( CUB.norPure{2}.LARGE ); % L -> large -> cor (veridical)
                elseif (idch == 10)
                    [norm_cube] = get_normCube( CUB.extPure{2}.LARGE );% XL -> large -> cor(veridical)
                % Episode type III: C = small and F = incorrect
                elseif (idch == 11)
                    [norm_cube] = get_normCube( CUB.extPure{3}.SMALL ); % XS -> small -> incor
                elseif (idch == 12)
                    [norm_cube] = get_normCube( CUB.norPure{3}.SMALL );% S -> small -> incor
                elseif (idch == 13)
                    [norm_cube] = get_normCube( CUB.extPure{7}.SMALL );% M -> small -> incor
                elseif (idch == 14) 
                    [norm_cube] = get_normCube( CUB.norPure{5}.SMALL ); % L -> small -> incor
                elseif (idch == 15)
                    [norm_cube] = get_normCube( CUB.extPure{5}.SMALL );% XL -> small -> incor

                % Episode type IV: C = large and F = incorrect
                elseif (idch == 16)
                    [norm_cube] = get_normCube( CUB.extPure{5}.LARGE );% XS -> large -> incor
                elseif (idch == 17)
                    [norm_cube] = get_normCube( CUB.norPure{5}.LARGE );% S -> large -> incor
                elseif (idch == 18)
                    [norm_cube] = get_normCube( CUB.extPure{7}.LARGE );% M -> large -> incor
                elseif (idch == 19)   
                    [norm_cube] = get_normCube( CUB.norPure{3}.LARGE );% L -> large -> incor
                elseif (idch == 20)
                    [norm_cube] = get_normCube( CUB.extPure{3}.LARGE );% XL -> large -> incor
                end
                
                episodefreq = NaN(1, nStim); 
                for iS = 1:nStim % XS to XL
                    iCh = 1; % small
                    iFb = 2; % correct
                    episodefreq(1, iS) = norm_cube(iS , iCh, iFb);
                end
                chunk1 = episodefreq; 

                for iS = 1:nStim % XS to XL
                    iCh = 2; % large
                    iFb = 2; % correct
                    episodefreq(1, iS) = norm_cube(iS , iCh, iFb);
                end
                chunk2 = episodefreq; 



                for iS = 1:nStim % XS to XL
                    iCh = 1; % small
                    iFb = 1; % incorrect
                    episodefreq(1, iS) = norm_cube(iS , iCh, iFb);
                end
                chunk3 = episodefreq; 


                for iS = 1:nStim % XS to XL
                    iCh = 2; % large
                    iFb = 1; % incorrect
                    episodefreq(1, iS) = norm_cube(iS , iCh, iFb);
                end
                chunk4 = episodefreq; 

            stack_idx = find( axis_id == axis_id(id) ); 
            indChunk(:,stack_idx ) = [chunk1, NaN,  chunk2, NaN,  chunk3,  NaN, chunk4];

        end

    end
   
    chunk_divid = 1:3; 
    
    
    figure(14+which_data);clf

    h_img = imagesc(indChunk);colorbar; caxis([0 0.25]); colormap(pink)
    set(h_img, 'AlphaData', ~isnan(indChunk)); axis square;axis off;

    if (which_data == 1)
        title('Humans')
        allHmap = indChunk; 
    elseif (which_data == 2)
        title('W-model')
        allWmap = indChunk; 
    elseif (which_data == 3)
        title('V-model')
        allVmap = indChunk; 
    end

end


figure(161);clf; % Fig S6

subplot(2,1,1); % Value-updating model (vs. Humans) 
    h_img = imagesc(allVmap - allHmap);colorbar;  colormap(jet)
    set(h_img, 'AlphaData', ~isnan(allHmap))
    axis square;axis off;
    caxis([-0.03 0.03]); title('V-model - Human', 'fontsize', 20); 

subplot(2,1,2); % World-updating model (vs. Humans) 
    h_img = imagesc(allWmap - allHmap);colorbar;  colormap(jet)
    set(h_img, 'AlphaData', ~isnan(allHmap))
    axis square;axis off;

    caxis([-0.03 0.03]);title('W-model - Human', 'fontsize', 20); 

%% %% %% %% rest
%% Get INDIVIDUALS  and statistical tests
clear indvOBS ALL_indvOBS indvWMod ALL_indvWMod indvVMod ALL_indvVMod
numSub = size(CUB.extPure{2}.SMALL, 4);

% define 1-D axis of a matrix to plot (toi)
axis_id = [1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20]; 

for which_data = 1:3


    if (which_data == 1)
        path_str.data       = [path_str.currpath, '/data']; 
        load([path_str.data, '/Episode_freq_OBS.mat'])
    elseif (which_data == 2)
        path_str.data       = [path_str.currpath, '/results_expost']; 
        load([path_str.data,'/BMBU/Episode_freq_Wmodel_verS445_FF.mat'])
    elseif (which_data == 3)
        path_str.data       = [path_str.currpath, '/results_expost']; 
        load([path_str.data, '/RLVU/Episode_freq_Vmodel_verS_softbias_plus_SSpercpt.mat'])
    end
    if (which_toi == -1)
        CUB = toiminusone_Pcub;
    elseif (which_toi == +1)
        CUB = toiplusone_Pcub;
    end


    for iSub = 1 : numSub
        indChunk = NaN(size(axis_id, 2), size(axis_id, 2)); % 20-by-20 matrix for colored maps in Fig 7E, 7F
            for id = 1 : length(axis_id)

                if (isnan(axis_id(id)))
                    indChunk(:, id) = NaN(size(axis_id, 2), 1);

                else
                    idch = axis_id(id); 
                        % Episode type I: C = small and F = correct
                        if (idch == 1)
                            [norm_cube] = get_normCube( CUB.extPure{2}.SMALL(:,:,:, iSub) );% XS -> small -> cor(veridical) (when this toi clamped)
                        elseif (idch == 2)
                            [norm_cube] = get_normCube( CUB.norPure{2}.SMALL(:,:,:, iSub) );% S -> small -> cor(veridical)
                        elseif (idch == 3)
                            [norm_cube] = get_normCube( CUB.extPure{6}.SMALL(:,:,:, iSub) );% M -> small -> cor
                        elseif (idch == 4) 
                            [norm_cube] = get_normCube( CUB.norPure{4}.SMALL(:,:,:, iSub) );% L -> small -> cor
                        elseif (idch == 5)
                            [norm_cube] = get_normCube( CUB.extPure{4}.SMALL(:,:,:, iSub) );% XL -> small -> cor

                % Episode type II: C = large and F = correct
                        elseif (idch == 6)
                            [norm_cube] = get_normCube( CUB.extPure{4}.LARGE(:,:,:, iSub) );% XS -> large -> cor
                        elseif (idch == 7)
                            [norm_cube] = get_normCube( CUB.norPure{4}.LARGE(:,:,:, iSub) ); % S -> large -> cor
                        elseif (idch == 8)
                            [norm_cube] = get_normCube( CUB.extPure{6}.LARGE(:,:,:, iSub) ); % M -> large -> cor
                        elseif (idch == 9)   
                            [norm_cube] = get_normCube( CUB.norPure{2}.LARGE(:,:,:, iSub) ); % L -> large -> cor (veridical)
                        elseif (idch == 10)
                            [norm_cube] = get_normCube( CUB.extPure{2}.LARGE(:,:,:, iSub) );% XL -> large -> cor(veridical)
                % Episode type III: C = small and F = incorrect

                        elseif (idch == 11)
                            [norm_cube] = get_normCube( CUB.extPure{3}.SMALL(:,:,:, iSub) ); % XS -> small -> incor
                        elseif (idch == 12)
                            [norm_cube] = get_normCube( CUB.norPure{3}.SMALL(:,:,:, iSub) );% S -> small -> incor
                        elseif (idch == 13)
                            [norm_cube] = get_normCube( CUB.extPure{7}.SMALL(:,:,:, iSub) );% M -> small -> incor
                        elseif (idch == 14) 
                            [norm_cube] = get_normCube( CUB.norPure{5}.SMALL(:,:,:, iSub) ); % L -> small -> incor
                        elseif (idch == 15)
                            [norm_cube] = get_normCube( CUB.extPure{5}.SMALL(:,:,:, iSub) );% XL -> small -> incor

                % Episode type IV: C = large and F = incorrect

                        elseif (idch == 16)
                            [norm_cube] = get_normCube( CUB.extPure{5}.LARGE(:,:,:, iSub) );% XS -> large -> incor
                        elseif (idch == 17)
                            [norm_cube] = get_normCube( CUB.norPure{5}.LARGE(:,:,:, iSub) );% S -> large -> incor
                        elseif (idch == 18)
                            [norm_cube] = get_normCube( CUB.extPure{7}.LARGE(:,:,:, iSub) );% M -> large -> incor
                        elseif (idch == 19)   
                            [norm_cube] = get_normCube( CUB.norPure{3}.LARGE(:,:,:, iSub) );% L -> large -> incor
                        elseif (idch == 20)
                            [norm_cube] = get_normCube( CUB.extPure{3}.LARGE(:,:,:, iSub) );% XL -> large -> incor



                        end
                        episodefreq = NaN(1, nStim); 
                                    clear episodefreq
                                    for iS = 1:5 % XS to XL
                                    iCh = 1; % small
                                    iFb = 2; % correct
                                    episodefreq(iS) = norm_cube(iS , iCh, iFb);
                                    end
                                    chunk1 = episodefreq; 
                                   
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
%                 title('Humans')

 


                 ALL_indvOBS(:,:,:,iSub) = indChunk; 
            elseif (which_data == 2)
%                 title('W-model')


                  ALL_indvWMod(:,:,:,iSub) = indChunk; 
            elseif (which_data == 3)
%                     title('V-model')  

                      ALL_indvVMod(:,:,:,iSub) = indChunk; 
            end


    % 
    end
end

%% Fig 7E & Fig 7F
vm_bonferroniN = 400; 
wm_bonferroniN = 400; 


vm_bonferroniN = 396; 
wm_bonferroniN = 395; 


axis_id = [1,2,3,4,5, NaN,6,7,8,9,10, NaN, 11,12,13,14,15, NaN, 16,17,18,19,20]; 
ALL_raw_h= 1000*ones(size(allWmap) ); ALL_raw_p = 1000*ones(size(allWmap) ); 
vALL_raw_h= 1000*ones(size(allWmap) ); vALL_raw_p = 1000*ones(size(allWmap) ); 

for id = 1 : size(axis_id, 2)
    if (isnan(axis_id(id)))
ALL_raw_h(:, id) = NaN(size(axis_id));
 ALL_raw_p(:, id) = NaN(size(axis_id));
 
 
 vALL_raw_h(:, id) = NaN(size(axis_id));
 vALL_raw_p(:, id) = NaN(size(axis_id));

    else
        j = axis_id(id); 
    
        
        
    for id2 = 1: size(axis_id, 2)
          if (isnan(axis_id(id2)))
              ALL_raw_h(id2,id) = NaN; 
              ALL_raw_p(id2, id) = NaN; 
              vALL_raw_h(id2,id) = NaN; 
              vALL_raw_p(id2, id) = NaN; 
          else
              i = axis_id(id2);
              [ALL_raw_h(id2,id), ALL_raw_p(id2, id) ] = ttest(ALL_indvOBS(i,j,1,:)  , ALL_indvWMod(i,j,1,:));
              if (ALL_raw_p(id2, id) < 0.05 / wm_bonferroniN )
                  ALL_raw_h(id2,id) = 1; 
              elseif (ALL_raw_p(id2, id) >= 0.05 / wm_bonferroniN )
                  
                    ALL_raw_h(id2,id) = 0; 
              end
             
            
                 [vALL_raw_h(id2,id), vALL_raw_p(id2, id) ] = ttest(ALL_indvOBS(i,j,1,:)  , ALL_indvVMod(i,j,1,:));
                 
                 
                 
                   if (vALL_raw_p(id2, id) < 0.05 / vm_bonferroniN )
                  vALL_raw_h(id2,id) = 1; 
              elseif (vALL_raw_p(id2, id) >= 0.05 / vm_bonferroniN )
                  
                    vALL_raw_h(id2,id) = 0; 
                   end
              
                   
          end
        
    end
    
    end
end
which_data = 2; 
figure(71+which_data);clf

h_img = imagesc(ALL_raw_h);colorbar; caxis([0 1]); colormap([0.5 0.5 0.5; 0 0 0 ])
set(h_img, 'AlphaData', ~isnan(ALL_raw_h))
axis square;axis off;

which_data = 3
figure(71+which_data);clf

h_img = imagesc(vALL_raw_h);colorbar; caxis([0 1]); colormap([0.5 0.5 0.5; 0 0 0 ])
set(h_img, 'AlphaData', ~isnan(vALL_raw_h))
axis square;axis off;


