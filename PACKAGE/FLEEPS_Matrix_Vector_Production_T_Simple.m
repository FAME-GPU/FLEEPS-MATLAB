function vec_y = FLEEPS_Matrix_Vector_Production_T_Simple(vec_x,Nx,Ny,Nz,N,D_k)

        % Reshape the input vector to do the inverse fast Fourier transform
        mat_x = reshape(vec_x,Nx,Ny,Nz);
        % Inverse Fast Fourier Transform 
        vec_y = reshape(ifftn(mat_x),N,1);
        %
        vec_y = D_k.*vec_y;
        % Assemble output vector
        vec_y = sqrt(N)*vec_y;
end