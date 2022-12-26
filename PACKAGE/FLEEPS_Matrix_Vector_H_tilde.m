function vec_y =  FLEEPS_Matrix_Vector_H_tilde( vec_x, fnchand_R_sqrt_QQ0, fnchand_QQ0s_R_sqrt, fnchand_H_tilde_11, fnchand_H_tilde_12, fnchand_H_tilde_21, fnchand_H_tilde_22)


   vec_y  = fnchand_QQ0s_R_sqrt(vec_x);

   vec_y1 =  fnchand_H_tilde_11 ( vec_y(1 : end-3) ) + fnchand_H_tilde_12 ( vec_y (end-2 : end) );
   vec_y2 = fnchand_H_tilde_21 ( vec_y(1 : end-3) ) + fnchand_H_tilde_22 ( vec_y (end-2 : end) );
   vec_y = [vec_y1; vec_y2];

   vec_y = fnchand_R_sqrt_QQ0(vec_y);
end