function B = FLEEPS_Matrix_B_Isotropic_phononic( Par_mesh, Par_material)

N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.B.org_idx);
N_lambda    = length(Par_material.lambda_in);
% N_mu        = length(Par_material.mu_in);

if (N_material_handle~=N_lambda)
    fprintf('\n');
    warning('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
    lambda_in = Par_material.lambda_in(1)*ones(N_material_handle,1);
    mu_in = Par_material.mu_in(1)*ones(N_material_handle,1);
    mu_in_face_center = Par_material.mu_in(1)*ones(N_material_handle,1);
else
    lambda_in = Par_material.lambda_in;
    mu_in = Par_material.mu_in;
    mu_in_face_center = Par_material.mu_in;
end

B_lambda_x = Par_material.lambda_out*ones(N,1);  
B_lambda_y = Par_material.lambda_out*ones(N,1);  
lambda_z = Par_material.lambda_out*ones(N,1);
B_mu_x  = Par_material.mu_out*ones(N,1);  
B_mu_y  = Par_material.mu_out*ones(N,1);  
B_mu_z  = Par_material.mu_out*ones(N,1);
B_mu_x_face_center  = Par_material.mu_out*ones(N,1);  
B_mu_y_face_center  = Par_material.mu_out*ones(N,1);  
B_mu_z_face_center  = Par_material.mu_out*ones(N,1);
for i = 1:N_material_handle
    B_inout_lambda_1 = zeros(N,1); 
    B_inout_lambda_2 = zeros(N,1); 
    B_inout_lambda_3 = zeros(N,1);
    B_inout_mu_1 = zeros(N,1); 
    B_inout_mu_2 = zeros(N,1); 
    B_inout_mu_3 = zeros(N,1);
    B_inout_mu_1_face_center = zeros(N,1); 
    B_inout_mu_2_face_center = zeros(N,1); 
    B_inout_mu_3_face_center = zeros(N,1);

    B_inout_lambda_1(Par_material.B.org_idx{1}) = 1; 
    B_inout_lambda_2(Par_material.B.org_idx{1}) = 1; 
    B_inout_lambda_3(Par_material.B.org_idx{1}) = 1;
    B_inout_mu_1(Par_material.B.org_idx{1}) = 1; 
    B_inout_mu_2(Par_material.B.org_idx{1}) = 1; 
    B_inout_mu_3(Par_material.B.org_idx{1}) = 1;
    B_inout_mu_1_face_center(Par_material.B.face_x_idx{1}) = 1; 
    B_inout_mu_2_face_center(Par_material.B.face_y_idx{1}) = 1; 
    B_inout_mu_3_face_center(Par_material.B.face_z_idx{1}) = 1;
    
    B.B_inout_ele(:,i) = [B_inout_lambda_1;B_inout_lambda_2;B_inout_lambda_3];
    B.B_inout_mag(:,i) = [B_inout_mu_1;B_inout_mu_2;B_inout_mu_3];
    B.B_inout_mag_face_center(:,i) = [B_inout_mu_1_face_center;B_inout_mu_2_face_center;B_inout_mu_3_face_center];
    
    B_lambda_x(Par_material.B.org_idx{i}) = lambda_in(i);
    B_lambda_y(Par_material.B.org_idx{i}) = lambda_in(i);
    lambda_z(Par_material.B.org_idx{i}) = lambda_in(i);
    B_mu_x(Par_material.B.org_idx{i})  = mu_in(i);
    B_mu_y(Par_material.B.org_idx{i})  = mu_in(i);
    B_mu_z(Par_material.B.org_idx{i})  = mu_in(i);
    B_mu_x_face_center(Par_material.B.face_x_idx{i})  = mu_in_face_center(i);
    B_mu_y_face_center(Par_material.B.face_y_idx{i})  = mu_in_face_center(i);
    B_mu_z_face_center(Par_material.B.face_z_idx{i})  = mu_in_face_center(i);
end
B.B_lambda    = [ B_lambda_x; B_lambda_y; lambda_z ];
B.B_mu     = [ B_mu_x ; B_mu_y ; B_mu_z  ];
B.B_mu_face_center     = [ B_mu_x_face_center ; B_mu_y_face_center ; B_mu_z_face_center  ];
B.invB_lambda = 1./B.B_lambda;
B.invB_mu  = 1./B.B_mu;
B.invB_mu_face_center  = 1./B.B_mu_face_center;
end
