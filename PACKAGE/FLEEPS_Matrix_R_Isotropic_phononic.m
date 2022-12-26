function R = FLEEPS_Matrix_R_Isotropic_phononic( Par_mesh, Par_material)

N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.B.edge_x_idx);
N_rho         = length(Par_material.rho_in);

if (N_material_handle~=N_rho)
    fprintf('\n');
    warning('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
    rho_in = Par_material.rho_in(1)*ones(N_material_handle,1);
else
    rho_in = Par_material.rho_in;
end

R_rho_x = Par_material.rho_out*ones(N,1);
R_rho_y = Par_material.rho_out*ones(N,1);  
R_rho_z = Par_material.rho_out*ones(N,1);

for i = 1:N_material_handle
    R_inout_edge_1 = zeros(N,1); 
    R_inout_edge_2 = zeros(N,1); 
    R_inout_edge_3 = zeros(N,1);
    R_inout_edge_1(Par_material.B.edge_x_idx{1}) = 1; 
    R_inout_edge_2(Par_material.B.edge_y_idx{1}) = 1; 
    R_inout_edge_3(Par_material.B.edge_z_idx{1}) = 1;
    R.R_inout_ele(:,i) = [R_inout_edge_1; R_inout_edge_2; R_inout_edge_3];
    
    R_rho_x(Par_material.B.edge_x_idx{i}) = rho_in(i);
    R_rho_y(Par_material.B.edge_y_idx{i}) = rho_in(i);
    R_rho_z(Par_material.B.edge_z_idx{i}) = rho_in(i);
end
R.R   = [ R_rho_x; R_rho_y; R_rho_z ];
R.invR = 1./R.R;
R.R_sqrt = sqrt(R.R);
R.invR_sqrt = 1./R.R_sqrt;
end
