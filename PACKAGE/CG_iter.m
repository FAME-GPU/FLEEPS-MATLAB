function [CG_cpu_time, CG_err, CG_iter] = CG_iter ( vec_x, fnchand_Qs_invS,  fnchand_S_Q, fnchand_PLP)

    tol = 1e-12;
   
    itmax  = 10000; 

    vec_y_initial= fnchand_Qs_invS (vec_x);

    time_start_pcg = tic;
    
    [vec_y,~,~,CG_iter] = pcg(fnchand_PLP, vec_y_initial, tol, itmax);

    CG_cpu_time = toc(time_start_pcg);

    CG_err =  norm(fnchand_S_Q( fnchand_PLP(vec_y) - vec_y_initial ));
end