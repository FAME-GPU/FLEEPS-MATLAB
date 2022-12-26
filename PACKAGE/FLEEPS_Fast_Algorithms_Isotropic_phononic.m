function [ freq, check, check_eigs, cpu_time, LS_iter, LS_cpu_time,  Ls_err] = FLEEPS_Fast_Algorithms_Isotropic_phononic( grid_num, FLEEPS_mtx, FLEEPS_fnchand, eig_option )
    
    global inner_iter inner_cpu_time inner_count  inner_err 
    inner_iter = zeros(1, eig_option.itmax);
    inner_cpu_time = zeros(1, eig_option.itmax);
    inner_err        = zeros(1, eig_option.itmax);

    inner_count = 0;

    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    N  = Nx*Ny*Nz;

    eig_option.v0 = eval(eig_option.v0);
    
    cpu_time_start = tic;
    opt.isreal = 0;
    opt.issym  = 1;
%     opt.v0     = eig_option.v0;
%     opt.v0     = eig_option.v0 - .01*FLEEPS_fnchand.Y0A0Y0(eig_option.v0);
    
    if length( fieldnames(FLEEPS_mtx.Q)) == 4
        [ev,ew] = eigs( FLEEPS_fnchand.H_tilde, 3*N, eig_option.eigen_wanted, 'lm', opt );
    else
        [ev,ew] = eigs( FLEEPS_fnchand.invA, 3*N, eig_option.eigen_wanted, 'lm', opt );
    end

    cpu_time = toc(cpu_time_start);
    
    LS_iter     = inner_iter;
    LS_cpu_time = inner_cpu_time;
    Ls_err         = inner_err;

    clear('inner_iter','inner_cpu_time');

    ew = 1./diag(ew);
    [ew,idx] = sort(ew,'ascend');
    ev = ev(:,idx);

    for i = 1:eig_option.eigen_wanted
        check_eigs(i,1) = norm (FLEEPS_fnchand.invA(ev(:,i)) - 1./ew(i)*ev(:,i),inf) / norm(ew(i)*ev(:,i));
        check(i,1) = norm (FLEEPS_fnchand.A(ev(:,i)) - ew(i)*ev(:,i),inf) / norm(ew(i)*ev(:,i));
    end

   freq = sqrt(ew);
   
end

