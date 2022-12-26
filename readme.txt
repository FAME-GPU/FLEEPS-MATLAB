# FLEEPS-MATLAB
The goal of FLEEPS-MATLAB is to solve the energy band of phononic crystal and get the accurate information of its energy band as soon as possible. 

# Installation
Download FLEEPS-MATLAB from Github.

# Introduction
Information of linear system in this paper can be found in "DATA/".

# Operating Instructions:
(1)  First, open FLEEPS_Compute_BandStructure_XXX.m to set your parameters.  
The following is the description of parameters:  
1. Lattice settings  
Popt.mesh.grid_num: The grid numbers;  
Popt.lattice.lattice_type: default;  
Popt.lattice.settingformat{1}: default;  
Popt.material.periodic_judge: default;  
Popt.material.display_material: Set whether to output material information( on/off );  
Popt.material.display_grid: Set whether to output mesh information( on/off);  
Popt.material.data_name: Set the lattice material, you can choose any one in \PACKAGE \FLEEPS_Material \FLEEPS_Material_Locate_Parameter;  
Popt.material.sphere_radius: Set the radius of the ball;  
Popt.material.cylinder_radius: default;  
Popt.material.material_type: default;  
Popt.material.rho_in, Popt.material.rho_out: Set densities in different materials;  
Popt.material.C_l_in, Popt.material.C_l_out: Set pressure wave velocities in different materials;  
Popt.material.C_t_in, Popt.material.C_t_out: Set shear wave velocities in different materials;  
Popt.material.lambda_in, Popt.material.lambda_out, Popt.material.mu_in, Popt.material.mu_out: Set Lame constants in different materials.  
2. Other settings   
Popt.eig.eigen_wanted: Set the number of eigenvalues ​​you want;  
Popt.eig.v0: Set the initial vector for the Lanczos method;  
Popt.eig.m: Set the subspace size of the Lanczos method;  
Popt.eig.itmax: Set the maximum number of iterations;  
Popt.eig.tolerence: Set the precision of the eigensolver;  
Popt.linsys.v0: Set Initial values ​​for linear system solver;  
Popt.linsys.tolerence: Set the precision of the linear system solver;  
Popt.linsys.itmax: Set the maximum number of iterations for the linear system solver;  
Popt.recip_lattice.part_num: Set the number of divisions between adjacent nodes on the first Brillouin zone path;  
Popt.recip_lattice.path_string: default;  
Popt.boundary.bd_cond_x, Popt.boundary.bd_cond_y, Popt.boundary.bd_cond_z: default;  
Popt.FLEEPS_option.discrete_method: default;  
Popt.FLEEPS_option.core_type: default.  

(2) Then just run FLEEPS_Compute_BandStructure_XXX.m to compute the band structure.  

(3) Use FLEEPS_Plot_LS_Info_XXX.m to plot the information of linear system of FLEEPS.  
     The relevant codes for performance comparison between AMG and (W)SVD can be found in PACKAGES/FLEEPS_MtxFnchand_Generate.m (line 36 -- line 64).

# Support:
If you have any questions, please let us know, we are willing to help you.  
School of Mathematics, Southeast University, Nanjing, China.  
Xing-Long Lyu: lxl_math@seu.edu.cn  
Hao-Nan Yang: hn_yang@seu.edu.cn  