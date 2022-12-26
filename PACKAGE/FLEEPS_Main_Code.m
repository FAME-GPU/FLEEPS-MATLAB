function [ Freq_array, Info ] = FLEEPS_Main_Code( Par_mesh, Par_lattice, Par_material, Par_eig, Par_linsys, wave_vec_array, Par_FLEEPS_option )
%FLEEPS_Main_Code Compute the source-free Maxwell's equations of a photonic crystal.

warning off;

fprintf('--------------------------------Start FLEEPS program----------------------------------\n')
for i = 2 :  size(wave_vec_array,2)
    
    %% (Main Step 1) - Construct FLEEPS matrices and function handles
    wave_vec = wave_vec_array(:,i); % Choose wave vector in Brillouin zone path
    fprintf('Start FLEEPS(%s, wave vector = [%.2e,%.2e,%.2e])\n', Par_material.material_type, wave_vec)
    [ FLEEPS_mtx, FLEEPS_fnchand ] = FLEEPS_MtxFnchand_Generate( wave_vec, Par_mesh, Par_lattice, Par_material, Par_linsys, Par_FLEEPS_option);
    
    %% (Main Step 3) - Use fast algorithms to compute the smallest frequency for given wave vector, material type, and lattice type
    switch Par_FLEEPS_option.discrete_method
        case 'staggered_grid'
            [ Freq_array(:,i) , Info.check(:,i), Info.check_eigs(:,i), Info.cpu_time(i), Info.LS_iter{i}, Info.LS_cpu_time{i}, Info.LS_err{i}] = ...
                FLEEPS_Fast_Algorithms( Par_mesh, Par_material, Par_eig, FLEEPS_mtx, FLEEPS_fnchand, Par_FLEEPS_option.core_type);
            if norm( wave_vec) == 0
                num = size( Freq_array(:,i), 1 );
                Freq_array(:,i) = [ zeros(3,1) ; Freq_array(1:7 ,i)];
            end
    end
    fprintf('Done! Use time %.2f(sec.)\n',Info.cpu_time(i))
    fprintf('================================== %3.1f%% complete =====================================\n',100*(i/size(wave_vec_array,2)))
end

end


function [ Freq_array, check, check_eigs, cpu_time, LS_iter, LS_cpu_time, Ls_err] = ...
    FLEEPS_Fast_Algorithms( Par_mesh, Par_material, Par_eig, FLEEPS_mtx, FLEEPS_fnchand, core_type)
switch Par_material.material_type
    case 'isotropic'
        if strcmp(core_type, 'cpuArray') == 1
            [ Freq_array, check, check_eigs, cpu_time, LS_iter, LS_cpu_time,  Ls_err] = ...
                FLEEPS_Fast_Algorithms_Isotropic_phononic( Par_mesh.grid_num, FLEEPS_mtx, FLEEPS_fnchand, Par_eig );
        end
end
end