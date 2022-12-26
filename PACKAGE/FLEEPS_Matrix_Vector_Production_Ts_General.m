function vec_y = FLEEPS_Matrix_Vector_Production_Ts_General( ...
      vec_x, Nx, Ny, Nz, N,...
      G_kx, G_ky, invsqrt_N_G_kz,...
      idx_1,idx_2)
    global mpv_time_Ts
    time_start = tic;
    %% reshape
    vec_y = vec_x;
    vec_y = reshape(vec_y, Nx, Ny*Nz);
    %% FFT(X-direction)
    vec_y = G_kx .* vec_y;
    vec_y = fft( vec_y );
    %% reshape
    vec_y = vec_y(idx_1);
    vec_y = reshape(vec_y,Ny,Nx*Nz);
    %% FFT(Y-direction)
    vec_y = G_ky .* vec_y;
    vec_y = fft( vec_y );
    %% reshape
    vec_y = vec_y(idx_2);
    vec_y = reshape(vec_y, Nz, Nx*Ny);
    %% FFT(Z-direction)
    vec_y = invsqrt_N_G_kz .* vec_y;
    vec_y = fft( vec_y );
    %% Output
    vec_y = reshape(vec_y, N, 1);
    mpv_time_Ts(end+1) = toc(time_start);
end

% function vec_z = FAME_Matrix_Vector_Production_Ts_General(vec_x,Nx,Ny,Nz,N,D_kx,D_ky,D_kz)
%     vec_z = zeros(N,1);
% 
%     vec_x_temp = fft(sparse(1:Nx,1:Nx,conj(D_kx)) * reshape(vec_x(1:N), Nx, Nz*Ny));
% 
%     vec_y_temp = zeros(Ny, Nx*Nz);
%     for ii = 1:Nx
%         vec_y_temp(:,(ii-1)*Nz+1:ii*Nz) = sparse(1:Ny,1:Ny,conj(D_ky(:,ii))) * reshape(vec_x_temp(ii,:).', Ny, Nz);
%     end
% 
%     vec_y_temp = fft(vec_y_temp);
% 
%     vec_z_temp = zeros(Nz, Nx*Ny);
%     for ii = 1:Nx
%         vec_z_temp(:,(ii-1)*Ny+1:ii*Ny) = conj(D_kz(:,1:Ny,ii)).*(vec_y_temp(:,(ii-1)*Nz+1:ii*Nz).');
%     end
% 
%     vec_z(1:N) = reshape(fft(vec_z_temp),Nx*Ny*Nz,1) / sqrt(N);
% end