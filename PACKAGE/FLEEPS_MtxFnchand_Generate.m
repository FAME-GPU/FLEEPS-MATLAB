function [ FLEEPS_mtx, FLEEPS_fnchand ] = FLEEPS_MtxFnchand_Generate(wave_vec, Par_mesh, Par_lattice, Par_material, Par_linsys, FLEEPS_option)
    
    % Construct material-dependent matrices
    switch FLEEPS_option.discrete_method
        case 'staggered_grid'
            switch Par_material.material_type
                case 'isotropic'
                FLEEPS_mtx.B = FLEEPS_Matrix_B_Isotropic_phononic( Par_mesh, Par_material);  % for phononic crystal 
                FLEEPS_mtx.L = FLEEPS_Matrix_L(FLEEPS_mtx.B.B_lambda, FLEEPS_mtx.B.B_mu, FLEEPS_mtx.B.B_mu_face_center);  % extra: construct  L(lambda, mu)
                FLEEPS_mtx.R = FLEEPS_Matrix_R_Isotropic_phononic( Par_mesh, Par_material);
                fprintf('The Matrix R has been formed\n')  
            end
    end
 
    % Construct discrete partial differential operators
    [ FLEEPS_mtx.C, FLEEPS_mtx.Cs, FLEEPS_mtx.C_1, FLEEPS_mtx.C_2, FLEEPS_mtx.C_3 ] = ...
        FLEEPS_Matrix_Construct( wave_vec, Par_mesh.grid_num, Par_mesh.edge_len,Par_mesh.mesh_len, {'quasi_periodic','quasi_periodic','quasi_periodic'}, Par_lattice.lattice_type, Par_lattice.lattice_vec_a_comp, Par_lattice.lattice_constant_comp );
    [ FLEEPS_mtx.Lambdas ] = ...
        FLEEPS_Matrix_Lambdas( wave_vec, Par_mesh.grid_num, Par_mesh.mesh_len, Par_lattice.lattice_type, Par_lattice.lattice_constant_comp, Par_lattice.lattice_vec_a_comp );

    [FLEEPS_mtx.Q, FLEEPS_mtx.S, FLEEPS_mtx.P] = FLEEPS_Matrix_Q_S_P_simplify(wave_vec, Par_mesh.grid_num, FLEEPS_mtx.Lambdas, Par_material , FLEEPS_mtx.R, FLEEPS_mtx.C_1, FLEEPS_mtx.C_2, FLEEPS_mtx.C_3);
    
    fprintf('The weighted SVD of D has been completed\n')

    % Construct function handles
    FLEEPS_fnchand = [];
    Nx = Par_mesh.grid_num(1); Ny = Par_mesh.grid_num(2); Nz = Par_mesh.grid_num(3);
    N  = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);

    switch FLEEPS_option.core_type
        case 'cpuArray'
        FLEEPS_fnchand = FLEEPS_Fnchandle_Manager( FLEEPS_fnchand, Nx, Ny, Nz, N, Par_material.material_type, Par_lattice.lattice_type, FLEEPS_mtx.B, FLEEPS_mtx.L, FLEEPS_mtx.R, FLEEPS_mtx.Q, FLEEPS_mtx.S, FLEEPS_mtx.P, FLEEPS_mtx.Lambdas, Par_linsys);
    end
    
    
    %% test preconditioner for AMG and (W)SVD
    for k = 1: 2
        t_0 = cputime;
        [FLEEPS_mtx_test.Q, FLEEPS_mtx_test.S, FLEEPS_mtx_test.P] = FLEEPS_Matrix_Q_S_P_simplify(wave_vec, Par_mesh.grid_num, FLEEPS_mtx.Lambdas, Par_material , FLEEPS_mtx.R, FLEEPS_mtx.C_1, FLEEPS_mtx.C_2, FLEEPS_mtx.C_3);
        t_svd(k) = cputime - t_0;
    end
    time_svd = mean(t_svd);
    
    t_1 = cputime;
    % Construct function handles
    FLEEPS_fnchand_test = [];
    Nx = Par_mesh.grid_num(1); Ny = Par_mesh.grid_num(2); Nz = Par_mesh.grid_num(3);
    N  = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
    FLEEPS_fnchand_test = FLEEPS_Fnchandle_Manager( FLEEPS_fnchand_test, Nx, Ny, Nz, N, Par_material.material_type, Par_lattice.lattice_type, FLEEPS_mtx.B, FLEEPS_mtx.L, FLEEPS_mtx.R, FLEEPS_mtx_test.Q, FLEEPS_mtx_test.S, FLEEPS_mtx_test.P, FLEEPS_mtx.Lambdas, Par_linsys);
    
    ze = sparse(N, N);
    D = [blkdiag(FLEEPS_mtx.C_1, FLEEPS_mtx.C_2, FLEEPS_mtx.C_3), ...
        [ze, -FLEEPS_mtx.C_3', -FLEEPS_mtx.C_2'; ...
        -FLEEPS_mtx.C_3', ze, -FLEEPS_mtx.C_1'; ...
        -FLEEPS_mtx.C_2', -FLEEPS_mtx.C_1', ze] ];
    A =  D * FLEEPS_mtx.L * D' ;
    A = (A + A')/2;
    b = ones(3 * N, 1);
    [FLEEPS_mtx_test.CG_cpu_time, FLEEPS_mtx_test.CG_err, FLEEPS_mtx_test.CG_iter] = CG_iter( b, FLEEPS_fnchand_test.Qs_invS,  FLEEPS_fnchand_test.S_Q,  FLEEPS_fnchand_test.PLP);
    FLEEPS_mtx_test.CG_cpu_time = cputime - t_1;
    [FLEEPS_mtx_test.amg_CG_cpu_time, FLEEPS_mtx_test.multi_x, FLEEPS_mtx_test.multi_info] = amg(FLEEPS_fnchand_test.A_without_R, FLEEPS_fnchand_test.Qs_invS,A,b);

    LS_info = [FLEEPS_mtx_test.multi_info.itStep, FLEEPS_mtx_test.multi_info.solverTime, FLEEPS_mtx_test.amg_CG_cpu_time,  FLEEPS_mtx_test.CG_iter, FLEEPS_mtx_test.CG_cpu_time, time_svd];
    writematrix(LS_info, 'LS_info.txt');
end


function FLEEPS_fnchand = FLEEPS_Fnchandle_Manager( FLEEPS_fnchand, Nx, Ny, Nz, N, material_type, lattice_type, B, L, R, Q, S, P, Lambdas, Par_linsys)
   
    Par_linsys_v0   = eval(Par_linsys.v0);

    % material-dependent matrices
    switch material_type
        case {'isotropic','Isotropic'}
            FLEEPS_fnchand.B_lambda    = @(vec_x) B.B_lambda.*vec_x;
            FLEEPS_fnchand.B_mu     = @(vec_x) B.B_mu .*vec_x;
            
            FLEEPS_fnchand.B_iso    = @(vec_x) vec_x.*B.B_lambda;
            FLEEPS_fnchand.invB_iso = @(vec_x) vec_x.*B.invB_lambda;
    end

    % FFT/IFFT-based matrices
    switch lattice_type
        case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            FLEEPS_fnchand.T   = @(vec_x) FLEEPS_Matrix_Vector_Production_T_Simple ( vec_x, Nx, Ny, Nz, N, Lambdas.D_k );
            FLEEPS_fnchand.Ts  = @(vec_x) FLEEPS_Matrix_Vector_Production_Ts_Simple( vec_x, Nx, Ny, Nz, N, Lambdas.D_ks);
        otherwise
            FLEEPS_fnchand.T = @(vec_x) FLEEPS_Matrix_Vector_Production_T_General(  vec_x, Nx, Ny, Nz, N, sqrt(N)*Lambdas.F_kx, Lambdas.F_ky, Lambdas.F_kz, Lambdas.idx_col_in_T, Lambdas.idx_in_T_1, Lambdas.idx_in_T_2 );
            FLEEPS_fnchand.Ts  = @(vec_x) FLEEPS_Matrix_Vector_Production_Ts_General( vec_x, Nx, Ny, Nz, N, Lambdas.G_kx, Lambdas.G_ky, Lambdas.G_kz/sqrt(N), Lambdas.idx_in_Ts_1, Lambdas.idx_in_Ts_2);
    end
    FLEEPS_fnchand.T3  = @(vec_x) [ FLEEPS_fnchand.T(  vec_x(    1:  N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(  N+1:2*N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(2*N+1:3*N,1) )];
    FLEEPS_fnchand.T3s = @(vec_x) [ FLEEPS_fnchand.Ts( vec_x(    1:  N,1) ) ;
                                  FLEEPS_fnchand.Ts( vec_x(  N+1:2*N,1) ) ;
                                  FLEEPS_fnchand.Ts( vec_x(2*N+1:3*N,1) )];

    FLEEPS_fnchand.T6  = @(vec_x) [ FLEEPS_fnchand.T(  vec_x(    1:  N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(  N+1:2*N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(2*N+1:3*N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(3*N+1:4*N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(4*N+1:5*N,1) ) ;
                                  FLEEPS_fnchand.T(  vec_x(5*N+1:6*N,1) )];
    FLEEPS_fnchand.T6s = @(vec_x) [ FLEEPS_fnchand.Ts( vec_x(    1:  N,1) ) ;
                                  FLEEPS_fnchand.Ts( vec_x(  N+1:2*N,1) ) ;
                                  FLEEPS_fnchand.Ts( vec_x(2*N+1:3*N,1) ) ;
                                  FLEEPS_fnchand.Ts(  vec_x(3*N+1:4*N,1) ) ;
                                  FLEEPS_fnchand.Ts(  vec_x(4*N+1:5*N,1) ) ;
                                  FLEEPS_fnchand.Ts(  vec_x(5*N+1:6*N,1) )];
    % SVD matrices
    FLEEPS_fnchand.P   = @(vec_x) FLEEPS_fnchand.T3( Lambdas.Pi_P  * vec_x );
    FLEEPS_fnchand.Pr  = @(vec_x) FLEEPS_fnchand.T3( Lambdas.Pi_Pr * vec_x );
    FLEEPS_fnchand.Q   = @(vec_x) FLEEPS_fnchand.T3( Lambdas.Pi_Q  * vec_x );
    FLEEPS_fnchand.Qr  = @(vec_x) FLEEPS_fnchand.T3( Lambdas.Pi_Qr * vec_x );
    FLEEPS_fnchand.Ps  = @(vec_x) Lambdas.Pi_Ps  * FLEEPS_fnchand.T3s( vec_x );
    FLEEPS_fnchand.Prs = @(vec_x) Lambdas.Pi_Prs * FLEEPS_fnchand.T3s( vec_x );
    FLEEPS_fnchand.Qs  = @(vec_x) Lambdas.Pi_Qs  * FLEEPS_fnchand.T3s( vec_x );
    FLEEPS_fnchand.Qrs = @(vec_x) Lambdas.Pi_Qrs * FLEEPS_fnchand.T3s( vec_x );

    FLEEPS_fnchand.R         = @(vec_x) R.R.*vec_x;
    FLEEPS_fnchand.R_sqrt    = @(vec_x) R.R_sqrt.*vec_x;
    FLEEPS_fnchand.invR      = @(vec_x) R.invR.*vec_x;
    FLEEPS_fnchand.invR_sqrt = @(vec_x) R.invR_sqrt.*vec_x;

    FLEEPS_fnchand.L         = @(vec_x) L*vec_x;
    

    FLEEPS_fnchand.Q         = @(vec_x) FLEEPS_fnchand.T3(  Q.Q * vec_x );
    FLEEPS_fnchand.Qs        = @(vec_x) Q.Qs * FLEEPS_fnchand.T3s( vec_x );
    FLEEPS_fnchand.Q0        = @(vec_x) FLEEPS_fnchand.T3(  Q.Q0 * vec_x );
    FLEEPS_fnchand.Q0s       = @(vec_x) Q.Q0s * FLEEPS_fnchand.T3s( vec_x );
    FLEEPS_fnchand.QQ0       = @(vec_x) FLEEPS_fnchand.T3(  [Q.Q, Q.Q0] * vec_x );
    FLEEPS_fnchand.QQ0s      = @(vec_x) [Q.Qs; Q.Q0s] * FLEEPS_fnchand.T3s( vec_x );

    FLEEPS_fnchand.P         = @(vec_x) FLEEPS_fnchand.T6(  P.P * vec_x );
    FLEEPS_fnchand.Ps        = @(vec_x) P.Ps * FLEEPS_fnchand.T6s( vec_x );

    FLEEPS_fnchand.S         = @(vec_x) S.S * vec_x;
    FLEEPS_fnchand.invS      = @(vec_x) S.invS * vec_x;
    FLEEPS_fnchand.R_sqrt_Qs_invS = @(vec_x)    FLEEPS_fnchand.invS ( FLEEPS_fnchand.Qs ( FLEEPS_fnchand.R_sqrt ( vec_x ) ) );
    FLEEPS_fnchand.invS_Q_R_sqrt = @(vec_x)   FLEEPS_fnchand.R_sqrt  ( FLEEPS_fnchand.Q  ( FLEEPS_fnchand.invS ( vec_x ) ) );

    FLEEPS_fnchand.Qs_invS = @(vec_x)    FLEEPS_fnchand.invS ( FLEEPS_fnchand.Qs ( vec_x ) );
    FLEEPS_fnchand.invS_Q = @(vec_x)    FLEEPS_fnchand.Q  ( FLEEPS_fnchand.invS ( vec_x ) ) ;

    FLEEPS_fnchand.S_Q    = @(vec_x)   FLEEPS_fnchand.Q( FLEEPS_fnchand.S  (  vec_x  ) );


    FLEEPS_fnchand.Qs_S    = @(vec_x)   FLEEPS_fnchand.S( FLEEPS_fnchand.Qs  (  vec_x  ) );
    FLEEPS_fnchand.S_Q_invR    = @(vec_x) FLEEPS_fnchand.invR (FLEEPS_fnchand.Q  ( FLEEPS_fnchand.S  ( vec_x ) ) );
    FLEEPS_fnchand.S_Qs_invR_sqrt = @(vec_x) FLEEPS_fnchand.S( FLEEPS_fnchand.Qs( FLEEPS_fnchand.invR_sqrt(vec_x) ) ); 
    FLEEPS_fnchand.invR_sqrt_Q_S    = @(vec_x) FLEEPS_fnchand.invR_sqrt( FLEEPS_fnchand.Q( FLEEPS_fnchand.S(vec_x) ) ); 

    FLEEPS_fnchand.PLP       = @(vec_x) FLEEPS_fnchand.Ps( FLEEPS_fnchand.L  ( FLEEPS_fnchand.P(vec_x) ) );
    FLEEPS_fnchand.A         = @(vec_x) FLEEPS_fnchand.invR_sqrt_Q_S  ( FLEEPS_fnchand.PLP  ( FLEEPS_fnchand.S_Qs_invR_sqrt (vec_x) ) );
    FLEEPS_fnchand.invA  = @(vec_x) FLEEPS_Matrix_Vector_Production_Isotropic_invA( vec_x, FLEEPS_fnchand.R_sqrt_Qs_invS,  FLEEPS_fnchand.invS_Q_R_sqrt, FLEEPS_fnchand.PLP, Par_linsys.tolerence, Par_linsys.itmax, Par_linsys_v0);
    FLEEPS_fnchand.A_without_R = @(vec_x) FLEEPS_fnchand.S_Q ( FLEEPS_fnchand.PLP ( FLEEPS_fnchand.Qs_S (vec_x ) ) );


    
    FLEEPS_fnchand.inv_S_PLP_S = @(vec_x) FLEEPS_Matrix_Vector_inv_S_PLP_S( vec_x, FLEEPS_fnchand.invS, FLEEPS_fnchand.PLP, Par_linsys.tolerence, Par_linsys.itmax);
    FLEEPS_fnchand.Q0RQ0     = @(vec_x) FLEEPS_fnchand.Q0s( FLEEPS_fnchand.R ( FLEEPS_fnchand.Q0( vec_x ) ) );
    FLEEPS_fnchand.inv_Q0RQ0 = @(vec_x) FLEEPS_Matrix_Vector_inv_Q0RQ0( vec_x, FLEEPS_fnchand.Q0RQ0 ); 
    FLEEPS_fnchand.Qs_R_Q0   = @(vec_x)  FLEEPS_fnchand.Qs( FLEEPS_fnchand.R ( FLEEPS_fnchand.Q0( vec_x ) ) );
    FLEEPS_fnchand.Q0s_R_Q   = @(vec_x)  FLEEPS_fnchand.Q0s( FLEEPS_fnchand.R ( FLEEPS_fnchand.Q( vec_x ) ) );

    FLEEPS_fnchand.H_tilde_11 = FLEEPS_fnchand.inv_S_PLP_S;
    FLEEPS_fnchand.H_tilde_12 = @(vec_x) -FLEEPS_fnchand.inv_S_PLP_S( FLEEPS_fnchand.Qs_R_Q0( FLEEPS_fnchand.inv_Q0RQ0 ( vec_x ) ) );
    FLEEPS_fnchand.H_tilde_21 = @(vec_x) -FLEEPS_fnchand.inv_Q0RQ0( FLEEPS_fnchand.Q0s_R_Q( FLEEPS_fnchand.inv_S_PLP_S ( vec_x ) ) );

    a0 = 1e12;
    FLEEPS_fnchand.H_tilde_22 = @(vec_x)  FLEEPS_fnchand.inv_Q0RQ0( FLEEPS_fnchand.Q0s_R_Q( FLEEPS_fnchand.inv_S_PLP_S( FLEEPS_fnchand.Qs_R_Q0( FLEEPS_fnchand.inv_Q0RQ0 ( vec_x ) ) ) ) ) + 1/a0 * FLEEPS_fnchand.inv_Q0RQ0 ( vec_x );

    FLEEPS_fnchand.R_sqrt_QQ0  = @(vec_x) FLEEPS_fnchand.R_sqrt ( FLEEPS_fnchand.QQ0( vec_x  ) );
    FLEEPS_fnchand.QQ0s_R_sqrt = @(vec_x) FLEEPS_fnchand.QQ0s (FLEEPS_fnchand.R_sqrt ( vec_x ) ) ;

    FLEEPS_fnchand.H_tilde = @(vec_x)  FLEEPS_Matrix_Vector_H_tilde( vec_x, FLEEPS_fnchand.R_sqrt_QQ0, FLEEPS_fnchand.QQ0s_R_sqrt, FLEEPS_fnchand.H_tilde_11, FLEEPS_fnchand.H_tilde_12, FLEEPS_fnchand.H_tilde_21, FLEEPS_fnchand.H_tilde_22);


    
    % NFSEP/NFGEP matrices        
    switch material_type
        case {'isotropic','Isotropic'}
            FLEEPS_fnchand.QBQ    = @(vec_x) FLEEPS_fnchand.Qrs( FLEEPS_fnchand.invB_iso( FLEEPS_fnchand.Qr(vec_x) ) );
            FLEEPS_fnchand.Ar  = @(vec_x) FLEEPS_Matrix_Vector_Production_Isotropic_Ar(vec_x, Lambdas.Sigma_r, FLEEPS_fnchand.QBQ);
            FLEEPS_fnchand.invAr  = @(vec_x) FLEEPS_Matrix_Vector_Production_Isotropic_invAr(vec_x, Lambdas.invSigma_r, FLEEPS_fnchand.QBQ, Par_linsys.tolerence, Par_linsys.itmax, Par_linsys_v0);
    end
end

