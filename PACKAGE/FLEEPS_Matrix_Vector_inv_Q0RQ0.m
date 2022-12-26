function vec_y = FLEEPS_Matrix_Vector_inv_Q0RQ0( vec_x, fnchand_Q0RQ0 )

  Mtx_s = zeros(3);
   for i = 1:3
       e  = zeros(3,1);
       e(i) = 1;
       Mtx_s(:,i) = fnchand_Q0RQ0(e);
   end

   vec_y = Mtx_s\vec_x;
end