function vec_y = FLEEPS_Matrix_Vector_inv_S_PLP_S (vec_x, fnchand_invS, fnchand_PLP, tol, itmax)


    vec_y= fnchand_invS (vec_x);
    
    [vec_y, ~, ~, ~] = pcg(fnchand_PLP, vec_y, tol, itmax);
%     [ vec_y, ~, ~, inner_iter(inner_count) ] = pcg(fnchand_PLP, vec_y, tol, itmax);
%     [ vec_y ] = pcg(fnchand_QBQ, vec_y, tol, itmax,[],[], v0);
    % inner_iter(inner_count) = 0;

    vec_y = fnchand_invS ( vec_y );
end