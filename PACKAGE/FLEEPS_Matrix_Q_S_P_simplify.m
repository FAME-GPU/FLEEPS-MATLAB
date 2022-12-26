function [Q ,S, P] = FLEEPS_Matrix_Q_S_P_simplify(wave_vec, grid_nums, Lambdas, Par_material, R, C_1, C_2, C_3)

n = grid_nums(1)*grid_nums(2)*grid_nums(3);

 for i = 1 : n
     ind_6n(1+6*(i-1) : 6*i, 1) = (i: n: i+5*n)'; 
 end


 for i = 1:n
     ind_3n(1+3*(i-1) : 3*i, 1) =  [i; n+i; 2*n+i];
 end


Lambda_1 = Lambdas.Lambda_x;
Lambda_2 = Lambdas.Lambda_y;
Lambda_3 = Lambdas.Lambda_z;
Lambda_1_conj = Lambda_1';
Lambda_2_conj = Lambda_2';
Lambda_3_conj = Lambda_3';


lambda = Par_material.lambda_in;
mu     = Par_material.mu_in;

coef_mu = sparse(1:3*n,1:3*n,2*mu,6*n,6*n) + sparse((1:3*n)+3*n, (1:3*n)+3*n,mu,6*n,6*n);
index_x_mu = [kron(ones(3,1),(1:n)'); kron(ones(3,1),(1+n:2*n)'); kron(ones(3,1),(1+2*n:3*n)')];
index_y_mu = kron(ones(3,1),(1:3*n)');

coef_mtx = coef_mu + sparse(index_x_mu,index_y_mu, lambda * ones(9*n,1) ,6*n,6*n);
coef_mtx_permutation =  coef_mtx(ind_6n, ind_6n);

M = coef_mtx_permutation(1:6,1:6);

R = chol(M);


if  norm(wave_vec) == 0

    for i = 1:n
        block3 = diag([Lambda_1(i), Lambda_2(i), Lambda_3(i)]);  
        block4 = [0, -Lambda_3_conj(i), -Lambda_2_conj(i); -Lambda_3_conj(i),0,-Lambda_1_conj(i);-Lambda_2_conj(i), -Lambda_1_conj(i),0];
         A = [block3, block4];
         A_RS = A*R';
         [U,S_temp,V_tilde] = svd(A_RS,'econ');

         V = R\V_tilde;

         % record ind_i to squezze zero eig value
         if norm(S_temp) == 0
             ind_i = i;
         end

         % storage all U
         ind_x_mtx = [i,i,i;i+n,i+n,i+n;i+2*n,i+2*n,i+2*n];
         ind_x_U(:,i) = ind_x_mtx(:);
         ind_y_mtx = [i,i+n,i+2*n;i,i+n,i+2*n;i,i+n,i+2*n];
         ind_y_U(:,i) = ind_y_mtx(:);
         U_all(:,i) = U(:);

         % storage all S
         ind_S(:,i) = [i, i+n, i+2*n]';
         S_all(:,i) = diag(S_temp);


         % storage all V
         ind_x_mtx = kron(ones(1,3),(i:n:i+5*n)');
         ind_x_V(:,i) = ind_x_mtx(:);
         ind_y_mtx = kron(ones(6,1),(i:n:i+2*n));
         ind_y_V(:,i) = ind_y_mtx(:);
         V_all(:,i) = V(:);
         
    end
      % form Q
         Q.Q  = sparse(ind_x_U, ind_y_U, U_all, 3*n, 3*n);
         % keep Q0
         Q.Q0 = Q.Q(:, [ind_i, ind_i+n, ind_i+2*n] );
         Q.Q0s = Q.Q0';
         Q.Q( :, [ind_i, ind_i+n, ind_i+2*n] ) = [];
     % form S
         
         S.S = sparse(ind_S, ind_S, S_all, 3*n, 3*n);
         S.S([ind_i, ind_i+n, ind_i+2*n], :) = [];
         S.S(:,[ind_i, ind_i+n, ind_i+2*n]) = [];

         S.invS = sparse(ind_S, ind_S, 1./S_all, 3*n, 3*n);
         S.invS([ind_i, ind_i+n, ind_i+2*n], :) = [];
         S.invS(:,[ind_i, ind_i+n, ind_i+2*n]) = [];

     % form P
         P.P = sparse(ind_x_V, ind_y_V, V_all, 6*n, 3*n);
         P.P( :, [ind_i, ind_i+n, ind_i+2*n] ) = [];

    
    
else

    for i = 1:n
        block3 = diag([Lambda_1(i), Lambda_2(i), Lambda_3(i)]);  
        block4 = [0, -Lambda_3_conj(i), -Lambda_2_conj(i); -Lambda_3_conj(i),0,-Lambda_1_conj(i);-Lambda_2_conj(i), -Lambda_1_conj(i),0];
         A = [block3, block4];
         A_RS = A*R';
         [U,S_temp,V_tilde] = svd(A_RS,'econ');

         V = R\V_tilde;
    
         % storage all U

         ind_x_mtx = [i,i,i;i+n,i+n,i+n;i+2*n,i+2*n,i+2*n];
         ind_x_U(:,i) = ind_x_mtx(:);
         ind_y_mtx = [i,i+n,i+2*n;i,i+n,i+2*n;i,i+n,i+2*n];
         ind_y_U(:,i) = ind_y_mtx(:);
         U_all(:,i) = U(:);

         % storage all S
         ind_S(:,i) = [i, i+n, i+2*n]';
         S_all(:,i) = diag(S_temp);


         % storage all V
         ind_x_mtx = kron(ones(1,3),(i:n:i+5*n)');
         ind_x_V(:,i) = ind_x_mtx(:);
         ind_y_mtx = kron(ones(6,1),(i:n:i+2*n));
         ind_y_V(:,i) = ind_y_mtx(:);
         V_all(:,i) = V(:);

    end
     % form Q
         Q.Q = sparse(ind_x_U, ind_y_U, U_all, 3*n, 3*n);

     % form S
         S.S    = sparse(ind_S, ind_S, S_all, 3*n, 3*n);

         S.invS = sparse(ind_S, ind_S, 1./S_all, 3*n, 3*n);
         % deal with zero
%          S.invS( find(S.invS == inf) ) = 0;

     % form P
         P.P = sparse(ind_x_V, ind_y_V, V_all, 6*n, 3*n);

end

 P.Ps = P.P';
 Q.Qs = Q.Q';
 
%  error = check_svd_3by6mtx(Lambdas, Q.Q, S.S, P.P, grid_nums, C_1, C_2, C_3);

end

%check for form_svd_3by6mtx
function error = check_svd_3by6mtx(Lambdas, Q,S,P,grid_nums, C_1, C_2, C_3)


Nx = grid_nums(1);
Ny = grid_nums(2);
Nz = grid_nums(3);
n = Nx*Ny*Nz;

  FLEEPS_fnchand.T   = @(vec_x) FLEEPS_Matrix_Vector_Production_T_Simple ( vec_x, Nx, Ny, Nz, n, Lambdas.D_k );
  FLEEPS_fnchand.Ts  = @(vec_x) FLEEPS_Matrix_Vector_Production_Ts_Simple( vec_x, Nx, Ny, Nz, n, Lambdas.D_ks);
x_vec1 = rand(6*n,1);
for i = 1:6
xs(1+n*(i-1):i*n,1) = FLEEPS_fnchand.Ts(x_vec1(1+n*(i-1):i*n,1));
end
temp = Q*S*P'*xs;
x_1 = temp(1:n,1);
x_2 = temp(n+1:2*n,1);
x_3 = temp(2*n+1:3*n,1);
result1 = FLEEPS_fnchand.T(x_1);
result2 = FLEEPS_fnchand.T(x_2);
result3 = FLEEPS_fnchand.T(x_3);
result = [result1;result2;result3];


% orign_mtx
blk1 = blkdiag(C_1,C_2,C_3);
blk2 = [zeros(n),-C_3',-C_2';-C_3',zeros(n),-C_1';-C_2',-C_1',zeros(n)];


origin = [blk1,blk2];

test_vec = origin*x_vec1;

error =norm(test_vec-result,"inf");

end


