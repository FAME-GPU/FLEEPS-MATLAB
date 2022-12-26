function   L = FLEEPS_Matrix_L(lambda,mu,mu_face_center)

n_3 = length(lambda);
n = n_3/3;
L_mu = sparse(1:n_3,1:n_3,2*mu,2*n_3,2*n_3) + sparse([1:n_3]+n_3,[1:n_3]+n_3,mu_face_center,2*n_3,2*n_3);

index_x = [kron(ones(3,1),[1:n]'); kron(ones(3,1),[1+n:2*n]'); kron(ones(3,1),[1+2*n:3*n]')];
index_y = kron(ones(3,1),[1:3*n]');
lambda_3l = kron(ones(3,1),lambda);
% Homo_mtx = sparse(index_x,index_y,[vector_all.line11;vector_all.line12;vector_all.line13;vector_all.line21;vector_all.line22;vector_all.line23;vector_all.line31;vector_all.line32;vector_all.line33],3*n,3*n);
L = L_mu + sparse(index_x,index_y,lambda_3l,2*n_3,2*n_3);

% constructure eig_line12
% eig_line(:,1) = lambda(1:n) + 2*mu(1:n);
% eig_line(:,2) = lambda(1+n:2*n);
% eig_line(:,3) = lambda(1+2*n:3*n);
% eig_line(:,4) = lambda(1:n) ;
% eig_line(:,5) = lambda(1+n:2*n) + 2*mu(1+n:2*n);
% eig_line(:,6) = lambda(1+2*n:3*n);
% eig_line(:,7) = lambda(1:n) ;
% eig_line(:,8) = lambda(1+n:2*n);
% eig_line(:,9) = lambda(1+2*n:3*n) + 2*mu(1+2*n:3*n);
% eig_line(:,10) = mu_face_center(1:n);
% eig_line(:,11) = mu_face_center(1+n:2*n);
% eig_line(:,12) = mu_face_center(1+2*n:3*n);
