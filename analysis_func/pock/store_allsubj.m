function store_allsubj(str_savefile)

    for iSub = 1 : 30
        load([str_savefile, '/Subj', num2str(iSub)], 'bootdata')



        obs_psepost.small_corr(iSub,:) = nanmean(bootdata.sSam, 3); 
        obs_psepost.large_corr(iSub,:) = nanmean(bootdata.lSam, 3);
        obs_psepost.small_incor(iSub,:) = nanmean(bootdata.sDiff, 3);
        obs_psepost.large_incor(iSub,:) = nanmean(bootdata.lDiff, 3);




    end
    temp = split(str_savefile, 'tmp/'); 
    if (strcmp(temp{2}, 'forward'))
        save([str_savefile, '/behav_prosp'], 'obs_psepost');
    elseif (strcmp(temp{2}, 'backward'))
        save([str_savefile, '/behav_retro'], 'obs_psepost');
    end
end
