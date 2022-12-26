function vec_y = FLEEPS_Matrix_Vector_Production_T_General(...
      vec_x, Nx, Ny, Nz, N,...
      sqrt_N_F_kx, F_ky, F_kz,...
      idx_col,idx_1,idx_2)
    global mpv_time_T
    time_start = tic;
    %% reshape
    vec_y = vec_x;
    vec_y = reshape(vec_y, Nz, Nx*Ny);
    %% IFFT(Z-direction)
    vec_y = ifft( vec_y );
    vec_y = F_kz .* vec_y;
    %% reshape
    vec_y = vec_y(idx_1);
    vec_y = reshape(vec_y, Ny, Nx*Nz);
    %% IFFT(Y-direction)
    vec_y = ifft( vec_y );
    vec_y = F_ky .* vec_y(:,idx_col);
    %% reshape
    vec_y = vec_y(idx_2);
    vec_y = reshape(vec_y, Nx, Ny*Nz);
    %% IFFT(X-direction)
    vec_y = ifft( vec_y );
    vec_y = sqrt_N_F_kx .* vec_y;
    %% Output
    vec_y = reshape(vec_y, N, 1);
    mpv_time_T(end+1) = toc(time_start);
end

% function vec_z = FAME_Matrix_Vector_Production_T_General(vec_x, Nx, Ny, Nz, N, D_kx, F_ky, F_kz)
%     vec_z = zeros(N,1);
%     vec_x_temp = reshape(vec_x(1:N), Nz, Nx*Ny);
%     vec_x_temp = F_kz.*(ifft(vec_x_temp));
% 
%     vec_y_temp = zeros(Ny, Nx*Nz);
% 
%     for ii = 1:Nx
%         vec_y_temp(:,(ii-1)*Nz+1:ii*Nz) = vec_x_temp(:,(ii-1)*Ny+1:ii*Ny).';
%     end
%     vec_y_temp = ifft(vec_y_temp);
% 
%     idx_col = zeros(Nx*Nz,1);
%     for ii = 1:Nz
%         idx_col((ii-1)*Nx+1:ii*Nx,1) = (ii:Nz:(Nx-1)*Nz+ii)';
%     end
%     vec_y_temp = F_ky.*vec_y_temp(:,idx_col);
% 
%     vec_z_temp = zeros(Nx, Ny*Nz);
%     for ii = 1:Nz
%         vec_z_temp(:,(ii-1)*Ny+1:ii*Ny) = vec_y_temp(:,(ii-1)*Nx+1:ii*Nx).';
%     end
%     vec_z(1:N) = sqrt(Nx * Ny * Nz) * reshape(sparse(1:Nx,1:Nx,D_kx)*ifft(vec_z_temp),Nx*Ny*Nz,1);
% end



