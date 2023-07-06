%% FLEEPS example for No186_Hexagonal_ZnS

clear
clc
close all

%% Load folder path of FLEEPS
addpath(genpath('PACKAGE/'));

%% Mesh setting
Popt.mesh.grid_num = [6,6,6]; % The grid numbers
%% Lattice settings
Popt.lattice.lattice_type     = 'user_defined';
Popt.lattice.settingformat{1} = 'Material Data';
% Popt.lattice.settingformat{2} = 'Constants form';
%% Material settings
Popt.material.periodic_judge   = 'on';
Popt.material.display_material = 'off';
Popt.material.display_grid     = 'off';
Popt.material.data_name       = 'No186_Hexagonal_ZnS';
Popt.material.sphere_radius   = 0.2; 
Popt.material.cylinder_radius = 0;
Popt.material.material_type   = 'isotropic';
Popt.material.rho_in  = 11357;
Popt.material.rho_out = 1180;
% Choice set pressure wave velocities and shear wave velocities  or  Lame constants
Popt.material.C_l_in  = 2158;
Popt.material.C_l_out = 2540;
Popt.material.C_t_in  = 860;
Popt.material.C_t_out = 1160;
[Popt.material.lambda_in , Popt.material.mu_in] = FLEEPS_Speed2LameCoef(Popt.material.C_l_in, Popt.material.C_t_in, Popt.material.rho_in);
[Popt.material.lambda_out, Popt.material.mu_out] = FLEEPS_Speed2LameCoef(Popt.material.C_l_out, Popt.material.C_t_out, Popt.material.rho_out);
Popt.eig.eigen_wanted  = 10; % This is the number of frequency (or eigenvalue) for each wave vector we want to plot on the graph,
Popt.eig.v0            = 'rand(3*N,1)'; % Lanczos initial
Popt.eig.m             = 30;          % Lanczos subspace size
Popt.eig.itmax         = 100;
Popt.eig.tolerence     = 1e-12;
Popt.linsys.v0        = 'rand(3*N,1)';
Popt.linsys.tolerence = 1e-12;
Popt.linsys.itmax     = 10000;
Popt.recip_lattice.part_num = 10; % The partitions of each path which is connected by two wave vectors
Popt.recip_lattice.path_string = 'default';
Popt.boundary.bd_cond_x = 'quasi_periodic';
Popt.boundary.bd_cond_y = 'quasi_periodic';
Popt.boundary.bd_cond_z = 'quasi_periodic';
Popt.FLEEPS_option.discrete_method = 'staggered_grid';
Popt.FLEEPS_option.core_type = 'cpuArray';

%% Generate modified lattice vectors and lattice constants for computing
[ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig, Par.linsys, Par.FLEEPS_option ] = ...
    FLEEPS_Parameter_Generator(Popt);

%% Generate wave vector array by partition the path of Brillouin zone 
[ Par.recip_lattice ] = ...
    FLEEPS_Parameter_Brillouin_Zone_Path( Par.recip_lattice.part_num, Par.lattice, Par.recip_lattice);

%% Locating indices for the material inside
[ Par.material.B.edge_x_idx, Par.material.B.edge_y_idx, Par.material.B.edge_z_idx, ... 
        Par.material.B.face_x_idx, Par.material.B.face_y_idx, Par.material.B.face_z_idx, Par.material.B.org_idx ] = ...
    FLEEPS_Material_Locate_Index( Par.mesh, Par.lattice, Par.material);

%% Start FLEEPS
[ omega_array, comput_info ] = ...
    FLEEPS_Main_Code( Par.mesh, Par.lattice, Par.material, Par.eig, Par.linsys, Par.recip_lattice.wave_vec_array, Par.FLEEPS_option );

%% Plot band structure
BS_ax = axes(figure);
FLEEPS_Plot_Band_Structure( Par.recip_lattice.path_string, Par.recip_lattice.part_num, omega_array/(2*pi), BS_ax );
