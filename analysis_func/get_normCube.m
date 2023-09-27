function [norm_cube] = get_normCube( input_epi )

    clamp_episode = input_epi;
    cube_cond = sum(clamp_episode, 4);
    norm_cube  = cube_cond/sum(cube_cond(:));
  
end