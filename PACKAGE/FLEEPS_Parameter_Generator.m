function [ Par_mesh, Par_lattice, Par_recip_lattice, Par_material, Par_eig, Par_linsys, Par_FLEEPS_option ] = FLEEPS_Parameter_Generator( Popt )
% This routine return the computational informations about given input:
%   grid_nums: a 1x3 array which contains grid numbers in terms of x-,y-, and z-axes, these must be positive integers.
%   lattice type:
%   lattice_constant:
%   Par_material:

grid_num      = Popt.mesh.grid_num;
lattice_type  = Popt.lattice.lattice_type;
settingformat = Popt.lattice.settingformat;
Par_material  = Popt.material;

if strcmp(settingformat{1},'Material Data')
    settingformat{2} = 'Constants form';
end

if isfield(Popt.lattice, 'lattice_constant') == 1
    lattice_constant_orig = Popt.lattice.lattice_constant;
else
    lattice_constant_orig = [];
end


switch settingformat{1}
    case {'Material Data'}
        % read material infomations from 'Par_material.data_name'
        data         = feval(Par_material.data_name,Par_material.sphere_radius,Par_material.cylinder_radius);
        lattice_type = data.lattice_type;
        Par_material.material_handle   = @(x,y,z,a1,a2,a3) FLEEPS_Material_Locate_Handle_User_Defined( x, y, z, a1, a2, a3, Par_material.data_name,Par_material.sphere_radius,Par_material.cylinder_radius);
        if isstruct(lattice_constant_orig) == 0
            lattice_constant_orig = data.lattice_constant;
        end
        Par_material.data = data;
        if length(Par_material.sphere_radius) ~= Par_material.data.material_num
            Par_material.sphere_radius = Par_material.sphere_radius(1)*ones(1,Par_material.data.material_num);
        end
        if length(Par_material.cylinder_radius) ~= Par_material.data.material_num
            Par_material.cylinder_radius = Par_material.cylinder_radius(1)*ones(1,Par_material.data.material_num);
        end
    otherwise
        error('Invalid Popt format!');
end


switch lattice_type
    case 'simple_cubic'
        lattice_constant_orig = FLEEPS_Parameter_Lattice_Constants_Format(lattice_type,lattice_constant_orig);
        
        Permutation    = [1 2 3];
        invPermutation = [1 2 3];
        Omega = eye(3);
        
        lattice_vec_a1 = lattice_constant_orig.a * [1;0;0];
        lattice_vec_a2 = lattice_constant_orig.b * [0;1;0];
        lattice_vec_a3 = lattice_constant_orig.c * [0;0;1];
        lattice_vec_a  = [lattice_vec_a1,lattice_vec_a2,lattice_vec_a3];
        lattice_vec_a_orig = lattice_vec_a;
        
        edge_len(1) = lattice_constant_orig.a;
        edge_len(2) = lattice_constant_orig.b;
        edge_len(3) = lattice_constant_orig.c;
        
        mesh_len     = [edge_len(1)/grid_num(1) edge_len(2)/grid_num(2) edge_len(3)/grid_num(3)];
        
        lattice_vec_a_comp    = lattice_vec_a;
        
        lattice_constant_comp.m1 = 0;
        lattice_constant_comp.m2 = 0;
        lattice_constant_comp.m3 = 0;
        lattice_constant_comp.m4 = 0;
        lattice_constant_comp.t1 = zeros(4,1);
        lattice_constant_comp.t2 = zeros(4,1);
        lattice_constant_comp.t3 = zeros(4,1);
        lattice_constant_comp.t4 = zeros(4,1);
        lattice_constant_comp.theta_1 = 0;
        lattice_constant_comp.theta_2 = 0;
        lattice_constant_comp.theta_3 = 0;
        lattice_constant_comp.length_a1 = 0;
        lattice_constant_comp.length_a2 = 0;
        lattice_constant_comp.length_a3 = 0;
    otherwise
        switch settingformat{2}
            case {'Constants form'}
                lattice_constant_orig = FLEEPS_Parameter_Lattice_Constants_Format(lattice_type,lattice_constant_orig);
                [lattice_vec_a_comp,lattice_vec_a_orig,lattice_constant_comp.length_a1,lattice_constant_comp.length_a2,lattice_constant_comp.length_a3,lattice_constant_comp.theta_1,lattice_constant_comp.theta_2,lattice_constant_comp.theta_3,Permutation] =...
                    FLEEPS_Parameter_Lattice_Vector(lattice_type,lattice_constant_orig);
        end
        
        Omega = lattice_vec_a_comp/lattice_vec_a_orig(:,Permutation);
        I = eye(3); P = I(:,Permutation); invP = P';
        invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
        
        edge_len(1) = lattice_vec_a_comp(1,1);
        edge_len(2) = lattice_vec_a_comp(2,2);
        edge_len(3) = lattice_vec_a_comp(3,3);
        mesh_len    = [edge_len(1)/grid_num(1) edge_len(2)/grid_num(2) edge_len(3)/grid_num(3)];
        
        [lattice_constant_comp.m1,lattice_constant_comp.m2,lattice_constant_comp.m3,lattice_constant_comp.m4,lattice_constant_comp.t1,lattice_constant_comp.t2,lattice_constant_comp.t3,lattice_constant_comp.t4] = ...
            FLEEPS_Parameter_Boundary_Point(lattice_vec_a_comp,lattice_constant_comp,mesh_len,grid_num);
end

Par_mesh.grid_num = grid_num;
Par_mesh.edge_len = edge_len;
Par_mesh.mesh_len = mesh_len;

Par_lattice.settingformat         = settingformat;
Par_lattice.lattice_type          = lattice_type;
Par_lattice.lattice_constant_comp = lattice_constant_comp;
Par_lattice.lattice_constant_orig = lattice_constant_orig;
Par_lattice.lattice_vec_a_comp    = lattice_vec_a_comp;
Par_lattice.lattice_vec_a_orig    = lattice_vec_a_orig;
Par_lattice.Permutation           = Permutation;
Par_lattice.invPermutation        = invPermutation;

Par_recip_lattice.part_num = Popt.recip_lattice.part_num;
Par_recip_lattice.lattice_vec_b_comp = inv(lattice_vec_a_comp');
Par_recip_lattice.lattice_vec_b_orig = inv(lattice_vec_a_orig');

Par_lattice.Omega = Omega;

Par_eig    = Popt.eig;
Par_linsys = Popt.linsys;

Par_FLEEPS_option = Popt.FLEEPS_option;
end
