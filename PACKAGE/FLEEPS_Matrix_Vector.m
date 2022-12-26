function result = FLEEPS_Matrix_Vector(L, eig_line, wave_vec,grid_nums)

[Q,S,P,P_eig_line] = svd_3by6mtx(wave_vec,grid_nums);
n = grid_nums(1)*grid_nums(2)*grid_nums(3); 
Nx = grid_nums(1);
Ny = grid_nums(2);
Nz = grid_nums(3);

edge_lens = ones(1,3);
bd_conds = {'quasi_periodic','quasi_periodic','quasi_periodic'};
mesh_lens = 1./grid_nums;
lattice_type = 'simple_cubic';
lattice_constant = [];
pml_grid_numbers= [];
lattice_vec_a = speye(3);


Lambdas = FLEEPS_Matrix_Lambdas( wave_vec, grid_nums, mesh_lens, lattice_type, lattice_constant, lattice_vec_a);
FLEEPS_fnchand.T   = @(vec_x) FLEEPS_Matrix_Vector_Production_T_Simple ( vec_x, Nx, Ny, Nz, n, Lambdas.D_k );
FLEEPS_fnchand.Ts  = @(vec_x) FLEEPS_Matrix_Vector_Production_Ts_Simple( vec_x, Nx, Ny, Nz, n, Lambdas.D_ks);

% extra: matrix-vector multiply
x_vec = rand(3*n,1);

%step1 : 18line × vector
for i = 1 : 6
      temp(:,i) = P_eig_line(:,1+3*(i-1)).*x_vec(1:n) +  P_eig_line(:,2+3*(i-1)).*x_vec(1+n:2*n) + P_eig_line(:,3*i).*x_vec(1+2*n:3*n);
end 


%step2: T × vector
for i = 1 : 6
    temp(:,i) = FLEEPS_fnchand.T(temp(:,i));
end


%step3: L mul vector
for i = 1 : 6
    if i <4
      temp2(:,i) = eig_line(:,1+3*(i-1)).*temp(:,1) +  eig_line(:,2+3*(i-1)).*temp(:,2) + eig_line(:,3*i).*temp(:,3);
    else
      temp2(:,i) = eig_line(:,6+i).*temp(:,i);
    end
end 



%step4: Ts mul vector
for i = 1 : 6
    temp2(:,i) = FLEEPS_fnchand.Ts(temp2(:,i));
end


%step5: 
for i = 1 : 3
    temp_result(:,i) =  conj(P_eig_line(:,i)).*temp2(:,1) + conj(P_eig_line(:,3+i)).*temp2(:,2) + conj(P_eig_line(:,6+i)).*temp2(:,3) + conj(P_eig_line(:,9+i)).*temp2(:,4) + conj(P_eig_line(:,12+i)).*temp2(:,5) + conj(P_eig_line(:,15+i)).*temp2(:,6);
end 

result = temp_result(:);
end


