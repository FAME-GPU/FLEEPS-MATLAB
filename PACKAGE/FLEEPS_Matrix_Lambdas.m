function Lambdas = FLEEPS_Matrix_Lambdas( wave_vec, grid_nums, mesh_lens, lattice_type, lattice_constant, lattice_vec_a )

    Nx = grid_nums(1); Ny = grid_nums(2); Nz = grid_nums(3);
   
    theta_x = dot(wave_vec,lattice_vec_a(:,1));
    theta_y = dot(wave_vec,lattice_vec_a(:,2));
    theta_z = dot(wave_vec,lattice_vec_a(:,3));
    switch lattice_type
        case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [Lambdas.Lambda_x, Lambdas.Lambda_y, Lambdas.Lambda_z, Lambdas.Lambda_q, Lambdas.Sigma, Lambdas.Sigma_r, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.Pi_Q, Lambdas.Pi_Qr, Lambdas.Pi_P, Lambdas.Pi_Pr, Lambdas.D_k, Lambdas.D_ks] =  ...
                Construct_Lambdas_Simple( wave_vec, Nx, Ny, Nz, mesh_lens, theta_x, theta_y, theta_z );
        otherwise
            if lattice_vec_a(1,2) > 0
                c1 = (   - lattice_constant.m1 / Nx );
            elseif lattice_vec_a(1,2) < 0
                c1 = ( 1 - lattice_constant.m1 / Nx );
            end
            a2_hat = c1*lattice_vec_a(:,1) + lattice_vec_a(:,2);
            a3_hat =    lattice_vec_a(:,3) - (lattice_constant.m3/Ny)*a2_hat - (lattice_constant.m2/Nx)*lattice_vec_a(:,1) + lattice_constant.t1;
                
            theta_xy  = wave_vec'*(a2_hat);
            theta_xyz = wave_vec'*(a3_hat);
            [ Lambdas.Lambda_x, Lambdas.Lambda_y, Lambdas.Lambda_z, Lambdas.Lambda_q, Lambdas.Sigma, Lambdas.Sigma_r,...
              Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_kx, Lambdas.F_ky, Lambdas.F_kz, Lambdas.G_kx, Lambdas.G_ky, Lambdas.G_kz,...
              Lambdas.idx_col_in_T, Lambdas.idx_in_T_1, Lambdas.idx_in_T_2, Lambdas.idx_in_Ts_1, Lambdas.idx_in_Ts_2,...
              Lambdas.Pi_Q, Lambdas.Pi_Qr, Lambdas.Pi_P, Lambdas.Pi_Pr] =  ...
                Construct_Lambdas_General( wave_vec, Nx, Ny, Nz, mesh_lens, theta_x, theta_xy, theta_xyz, lattice_constant);
    end
    
    Lambdas.invSigma_r = 1./Lambdas.Sigma_r;
    
    Lambdas.Pi_Qrs = Lambdas.Pi_Qr';
    Lambdas.Pi_Prs = Lambdas.Pi_Pr';
    Lambdas.Pi_Qs  = Lambdas.Pi_Q';
    Lambdas.Pi_Ps  = Lambdas.Pi_P';
    
    if norm(wave_vec) == 0
        Lambdas.Lambda_x = [0; Lambdas.Lambda_x ];
        Lambdas.Lambda_y = [0; Lambdas.Lambda_y ];
        Lambdas.Lambda_z = [0; Lambdas.Lambda_z ];
    end
    e = ones(Nx*Ny*Nz,1);
    Lambdas.Lambda_fem = (1/6)*[(mesh_lens(1)*Lambdas.Lambda_x + 2*e).*conj(mesh_lens(1)*Lambdas.Lambda_x + 2*e) + 2*e;
                                (mesh_lens(2)*Lambdas.Lambda_y + 2*e).*conj(mesh_lens(2)*Lambdas.Lambda_y + 2*e) + 2*e;
                                (mesh_lens(3)*Lambdas.Lambda_z + 2*e).*conj(mesh_lens(3)*Lambdas.Lambda_z + 2*e) + 2*e];
end

function [Lambda_x, Lambda_y, Lambda_z, Lambda_q, Sigma, Sigma_r, D_kx, D_ky, D_kz, Pi_Q, Pi_Qr, Pi_P, Pi_Pr, D_k, D_ks ] = ...
            Construct_Lambdas_Simple( wave_vec, Nx, Ny, Nz, mesh_lens, theta_x, theta_y, theta_z)
    e_x = ones(Nx,1);
    e_y = ones(Ny,1);
    e_z = ones(Nz,1);
    mx   = (0:Nx-1)';
    my   = (0:Ny-1)';
    mz   = (0:Nz-1)';
    
    temp1    = exp(1i*2*pi*(mx+theta_x)/Nx);
    temp2    = exp(1i*2*pi*(my+theta_y)/Ny); 
    temp3    = exp(1i*2*pi*(mz+theta_z)/Nz); 
    Lambda1   = (temp1-1)/mesh_lens(1);
    Lambda2   = (temp2-1)/mesh_lens(2);
    Lambda3   = (temp3-1)/mesh_lens(3);

    Lambda_x = kron(  e_z,kron(e_y,Lambda1));
    Lambda_y = kron(  e_z,kron(Lambda2,e_x));
    Lambda_z = kron(  Lambda3,kron(e_y,e_x));
       %% Diagonal matries associated with Bloch assumption
    D_kx = exp(theta_x*2*pi*1i*mx/Nx );
    D_ky = exp(theta_y*2*pi*1i*my/Ny );
    D_kz = exp(theta_z*2*pi*1i*mz/Nz );
        %% Construct Pr and Qr
    alpha = 3.1;  beta = 5.1;   
    N = Nz*Ny*Nx;
    Nd = N;
    diag_position = 0;
    if norm(wave_vec) == 0
        Nd = Nd-1;
        Lambda_x  = Lambda_x(2:end);
        Lambda_y  = Lambda_y(2:end);
        Lambda_z  = Lambda_z(2:end);
        diag_position = -1;
    end    
    Lambda_q = conj(Lambda_x).*Lambda_x + conj(Lambda_y).*Lambda_y + conj(Lambda_z).*Lambda_z;
    Lambda_p = [ beta*Lambda_z -       Lambda_y; 
                      Lambda_x - alpha*Lambda_z; 
                alpha*Lambda_y -  beta*Lambda_x];
    Lambda_q_sqrt = sqrt(Lambda_q);
    
    temp = conj(Lambda_p(1:Nd)).*Lambda_p(1:Nd) + conj(Lambda_p(Nd+1:2*Nd)).*Lambda_p(Nd+1:2*Nd) + conj(Lambda_p(2*Nd+1:3*Nd)).*Lambda_p(2*Nd+1:3*Nd);
    tempPi = temp.*Lambda_q;
    
    Pi_0 = [ spdiags( Lambda_x./Lambda_q_sqrt, diag_position,N, Nd ) ;
             spdiags( Lambda_y./Lambda_q_sqrt, diag_position,N, Nd ) ;
             spdiags( Lambda_z./Lambda_q_sqrt, diag_position,N, Nd ) ];
    D_k  = kron(D_kz,kron(D_ky,D_kx));
    D_ks = conj(D_k);
    if diag_position == -1
        temp_fft_ones = FLEEPS_Matrix_Vector_Production_Ts_Simple(ones(N,1)/sqrt(N),Nx,Ny,Nz,N,D_ks);
        Pi_0 = [Pi_0, kron(speye(3), temp_fft_ones ) ];
    end         
    Pi_1 = [ spdiags( ( (alpha*Lambda_y -  beta*Lambda_x).*conj(Lambda_y) - (      Lambda_x - alpha*Lambda_z).*conj(Lambda_z) ) ./ sqrt(tempPi), diag_position,N, Nd );
             spdiags( ( ( beta*Lambda_z -       Lambda_y).*conj(Lambda_z) - (alpha*Lambda_y -  beta*Lambda_x).*conj(Lambda_x) ) ./ sqrt(tempPi), diag_position,N, Nd );
             spdiags( ( (      Lambda_x - alpha*Lambda_z).*conj(Lambda_x) - ( beta*Lambda_z -       Lambda_y).*conj(Lambda_y) ) ./ sqrt(tempPi), diag_position,N, Nd )];

    Pi_2 = [ spdiags( ( beta*conj(Lambda_z) -       conj(Lambda_y))./sqrt(temp), diag_position,N, Nd ) ;
             spdiags( (      conj(Lambda_x) - alpha*conj(Lambda_z))./sqrt(temp), diag_position,N, Nd ) ;
             spdiags( (alpha*conj(Lambda_y) -  beta*conj(Lambda_x))./sqrt(temp), diag_position,N, Nd ) ];


    Sigma_r = [Lambda_q_sqrt;Lambda_q_sqrt];
    Pi_Qr   = sparse([      Pi_1 ,       Pi_2]);
    Pi_Pr   = sparse([-conj(Pi_2), conj(Pi_1)]);
    Sigma   = [zeros(3*N-2*Nd,1);Lambda_q_sqrt;Lambda_q_sqrt];
    Pi_Q    = sparse([       Pi_0,      Pi_1 ,       Pi_2]);
    Pi_P    = sparse([ conj(Pi_0),-conj(Pi_2), conj(Pi_1), ]);
end

function [Lambda_x, Lambda_y, Lambda_z, Lambda_q, Sigma, Sigma_r,...
          D_kx, D_ky, D_kz, F_kx, F_ky, F_kz, G_kx, G_ky, G_kz,...
          idx_col_in_T, idx_in_T_1, idx_in_T_2, idx_in_Ts_1, idx_in_Ts_2,...
          Pi_Q, Pi_Qr, Pi_P, Pi_Pr ] = ...
            Construct_Lambdas_General(wave_vec, Nx, Ny, Nz, mesh_lens, theta_x, theta_xy, theta_xyz, lattice_constant)
    m1 = lattice_constant.m1;  
    m2 = lattice_constant.m2;
    m3 = lattice_constant.m3;
    
    N = Nx*Ny*Nz;
    
    n1 = Nx;
    n2 = Ny;
   
    e_x  = ones(Nx,1);
    e_y  = ones(Ny,1);
    e_z  = ones(Nz,1);
    mx   = (0:Nx-1)';
    my   = (0:Ny-1)';
    mz   = (0:Nz-1)';
    mxy  = kron((-m1/Nx)*mx,e_y) + kron(e_x,my); 
    mxyz = ( (m1*m3)/(n1*n2) - ( m2/n1) ) * kron(  mx, kron( e_z, e_y) ) + ...
                                 (-m3/n2) * kron( e_x, kron(  my, e_z) ) + ...
                                            kron( e_x, kron( e_y,  mz) );

    temp1     = exp(1i*2*pi*(mx   + theta_x  )/Nx);
    temp12    = exp(1i*2*pi*(mxy  + theta_xy )/Ny); 
    temp123   = exp(1i*2*pi*(mxyz + theta_xyz)/Nz); 
    Lambda1   = (temp1   - 1)/mesh_lens(1);
    Lambda12  = (temp12  - 1)/mesh_lens(2);
    Lambda123 = (temp123 - 1)/mesh_lens(3);

    Lambda_x = kron( Lambda1 , kron(e_y,e_z));
    Lambda_y = kron( Lambda12, e_z);
    Lambda_z = Lambda123;
    %% Diagonal matries associated with Bloch assumption
    D_kx = exp(theta_x*2*pi*1i*mx/ Nx );
    D_ky = zeros(Ny,Nx);
    for jj = 1:Nx
        phi_jx     = 1i*2*pi * (theta_xy + (jj-1)*(-m1/Nx) ) / Ny;
        D_ky(:,jj) = exp( phi_jx * my );
    end
    D_kz = zeros(Nz,Ny,Nx);
    for jj = 1:Nx
        for ll = 1:Ny
            psi_k         = 1i*2*pi * ( theta_xyz + (ll-1)*(-m3/n2) + (jj-1)*((m1*m3)/(n1*n2) - (m2/n1)) ) / Nz; 
            D_kz(:,ll,jj) = exp( psi_k * mz );
        end
    end
    
    %% matrices used in T and Ts
    ones_Nx  = ones(1,Nx   ); ones_Ny  = ones(1,Ny   ); ones_Nz = ones(1,Nz);
    ones_Nxy = ones(1,Nx*Ny); ones_Nyz = ones(1,Ny*Nz); 

    idx_col_in_T = kron( ones_Nz, (0:Nx-1)*Nz ) + kron( 1:Nz, ones_Nx );
    
    tmp = kron(1:Nz,ones_Ny) + kron(ones_Nz,(0:Ny-1)*Nz);
    idx_in_T_1  = kron((0:Nx-1)*Ny*Nz,ones_Nyz) + kron( ones_Nx, tmp);
    tmp = kron(1:Ny,ones_Nx) + kron(ones_Ny,(0:Nx-1)*Ny);
    idx_in_T_2  = kron((0:Nz-1)*Nx*Ny,ones_Nxy) + kron( ones_Nz, tmp);

    idx_in_Ts_1 = kron(ones_Nx,(0:Ny*Nz-1)*Nx) + kron(1:Nx,ones_Nyz);
    tmp = kron(1:Ny,ones_Nz) + kron(ones_Ny,(0:Nz-1)*Ny);
    idx_in_Ts_2 = kron((0:Nx-1)*Ny*Nz,ones_Nyz) + kron( ones_Nx, tmp);

    F_kx  = kron(      D_kx , ones(1,Ny*Nz) );
    G_kx  =         kron( conj(D_kx), ones(1,Ny*Nz) );

    F_kz     = zeros(Nz, Nx*Ny);
    for ii = 1:Nx
        F_kz(:,(ii-1)*Ny+1:ii*Ny) = D_kz(:,:,ii);
    end
    G_kz = conj(F_kz);

    F_ky = kron(   ones(1,Nz)        , D_ky );
    G_ky = kron(   conj(D_ky), ones(1,Nz)   );
    %% Construct Pr and Qr
    diag_position = 0;
    Lambda_q = conj(Lambda_x).*Lambda_x + conj(Lambda_y).*Lambda_y + conj(Lambda_z).*Lambda_z;
    Lambda_s = Lambda_x + Lambda_y + Lambda_z;
    Lambda_p = conj(Lambda_s).*Lambda_s;
    if norm(wave_vec) < 1e-12
        Nd = N-1;
        Lambda_x  = Lambda_x(2:end);
        Lambda_y  = Lambda_y(2:end);
        Lambda_z  = Lambda_z(2:end);
        Lambda_q = Lambda_q(2:end);
        Lambda_s = Lambda_s(2:end);
        Lambda_p = Lambda_p(2:end);
        diag_position = -1;
    else
        Nd = N;
    end
    norm_Q1 = sqrt(3*Lambda_q.^2 - Lambda_q.*Lambda_p);
    norm_Q2 = sqrt(3*Lambda_q - Lambda_p);

    Pi_0 = [ spdiags( Lambda_x./sqrt(Lambda_q), diag_position, N, Nd ) ;
             spdiags( Lambda_y./sqrt(Lambda_q), diag_position, N, Nd ) ;
             spdiags( Lambda_z./sqrt(Lambda_q), diag_position, N, Nd )  ];
    if diag_position == -1
        temp_fft_ones = FLEEPS_Matrix_Vector_Production_Ts_General(ones(N,1)/sqrt(N),Nx,Ny,Nz,N,G_kx,G_ky,G_kz,idx_in_Ts_1,idx_in_Ts_2);
%         temp_fft_ones = FLEEPS_Matrix_Vector_Production_Ts_General(ones(N,1)/sqrt(N),Nx,Ny,Nz,N,D_kx,D_ky,D_kz);
        Pi_0 = [Pi_0, kron(speye(3), temp_fft_ones ) ];
    end         
    Pi_1 = [ spdiags( (Lambda_q - Lambda_x.*conj(Lambda_s))./norm_Q1, diag_position, N, Nd ) ;
             spdiags( (Lambda_q - Lambda_y.*conj(Lambda_s))./norm_Q1, diag_position, N, Nd ) ;
             spdiags( (Lambda_q - Lambda_z.*conj(Lambda_s))./norm_Q1, diag_position, N, Nd )  ];
    Pi_2 = [ spdiags( (conj(Lambda_z) - conj(Lambda_y))./norm_Q2, diag_position, N, Nd ) ;
             spdiags( (conj(Lambda_x) - conj(Lambda_z))./norm_Q2, diag_position, N, Nd ) ;
             spdiags( (conj(Lambda_y) - conj(Lambda_x))./norm_Q2, diag_position, N, Nd )  ];

    Sigma_r = [sqrt(Lambda_q);sqrt(Lambda_q)];
    Pi_Qr   = sparse([      Pi_1 ,       Pi_2]);
    Pi_Pr   = sparse([-conj(Pi_2), conj(Pi_1)]);
    Sigma   = [zeros(3*N-2*Nd,1);sqrt(Lambda_q);sqrt(Lambda_q)];
    Pi_Q    = sparse([       Pi_0,      Pi_1 ,       Pi_2]);
    Pi_P    = sparse([ conj(Pi_0),-conj(Pi_2), conj(Pi_1)]);
end

