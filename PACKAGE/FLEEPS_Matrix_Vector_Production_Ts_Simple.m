function vec_y = FLEEPS_Matrix_Vector_Production_Ts_Simple(vec_x,Nx,Ny,Nz,N,D_ks)

        vec_x = D_ks.*vec_x;
        % Reshape the input vector to do the fast Fourier transform
        mat_x = reshape(vec_x,Nx,Ny,Nz);
        % Fast Fourier Transform
        vec_y = reshape(fftn(mat_x),N,1);
        % Assemble output vector
        vec_y = (1/sqrt(N))*vec_y;   
end