function [ edge_x_idx, edge_y_idx, edge_z_idx, face_x_idx, face_y_idx, face_z_idx, org_idx ] = FLEEPS_Material_Locate_Index( Par_mesh, Par_lattice, Par_material )
% This programe return the indices corresponding to discrete points inside the material.
% Input: 
%    'grid_nums': must be a array which contains only 3 positive integers,
%           which stand for the discrete size for x-, y- and z- axis.
%    'lattice_type': must be a string array. It can be choosen as
%           'Simple_cubic', 'Face_centered_cubic' or 'Body_centered_cubic'
%    'lattice_constant': must be a positive real number.
%    'material_handle': it must be a function handle with the following
%           input and ouput forms:
%                index_in_material = material_handle(x,y,z,a1,a2,a3);
%           x,y,z are consider as the component of discrete points over the
%           parallel hexahedron spaned by a1, a2 and a3. For an example,
%           we may using following function handle to type spheres at
%           each corner of the parallel hexahedron corresponding to given a1, a2, a3.
%                    r = 0.4;
%                    sphere_handle = @(x,y,z,a1,a2,a3) find( (x).^2 + (y).^2 + (z).^2 < r^2 );
%           'material_handle' could be set as cell type parameter for two or more material. For example:
%                   gyroid_handle{1} = @(x,y,z,a1,a2,a3) FLEEPS_Material_Locate_Handle_Gyroid(  x,  y,  z, a1, a2, a3, 1.1, Popt.domain.lattice_constant );
%                   gyroid_handle{2} = @(x,y,z,a1,a2,a3) FLEEPS_Material_Locate_Handle_Gyroid( -x, -y, -z, a1, a2, a3, 1.1, Popt.domain.lattice_constant );
%    'display_material': default seting is 'off'. This function display the
%           material grid in a figure(12) if 'display_material' is put as a string
%           'on'.
% Output: the inside-material indices for Yee's grid, corresponding to these input.

    [ edge_center_x, edge_center_y, edge_center_z, face_center_x, face_center_y, face_center_z, origin_point_set] = ...
        FLEEPS_Matrix_Grid( Par_mesh.grid_num, Par_mesh.mesh_len );                              

    if isfield(Par_material,'color_map') == 0
        Par_material.color_map = { [15 131 225]/255, [200,40,45]/255, [238,126,2]/255, [75,165,102]/255, rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
    end
    
      %% Constructing the discretized grids on edge and face
    [ edge_x_idx ] = FLEEPS_Material_Locate( edge_center_x, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );
    [ edge_y_idx ] = FLEEPS_Material_Locate( edge_center_y, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );
    [ edge_z_idx ] = FLEEPS_Material_Locate( edge_center_z, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );
    
    [ face_x_idx ] = FLEEPS_Material_Locate( face_center_x, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );
    [ face_y_idx ] = FLEEPS_Material_Locate( face_center_y, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );
    [ face_z_idx ] = FLEEPS_Material_Locate( face_center_z, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );

    [  org_idx, reshaped_point_set  ] = FLEEPS_Material_Locate( origin_point_set, Par_lattice.lattice_vec_a_comp, Par_lattice.Omega, Par_lattice.Permutation, Par_material.material_handle, Par_material.periodic_judge );

    %% Plot Original material grids
    if iscell(org_idx) == 1
        material_num = length(org_idx);   
    else
        material_num   = 1;
        org_idx        = {org_idx};
        edge_x_idx      = {edge_x_idx};
        edge_y_idx      = {edge_y_idx};
        edge_z_idx      = {edge_z_idx};
        face_x_idx      = {face_x_idx};
        face_y_idx      = {face_y_idx};
        face_z_idx      = {face_z_idx};
    end
    
    if strcmp( Par_material.display_grid, 'on') || strcmp( Par_material.display_grid, 'On')
        hax_orig = axes(figure);
        cla(hax_orig);  hold(hax_orig,'on')
        for i = 1:material_num
            plot_material(hax_orig,reshaped_point_set, org_idx{i}, Par_lattice.Omega'*Par_lattice.lattice_vec_a_comp, Par_lattice.lattice_constant_comp, Par_material.color_map{i}, 'original');
        end
%         title('Material grids in original domain')
    end
    %% Plot Computational material grids 
    if strcmp( Par_material.display_grid, 'on') || strcmp( Par_material.display_grid, 'On')
        hax_comp = axes(figure);
        cla(hax_comp);  hold(hax_comp,'on');
        plotmesh3D(hax_comp,[0,Par_lattice.lattice_vec_a_comp(1,1)],...
                            [0,Par_lattice.lattice_vec_a_comp(2,2)],...
                            [0,Par_lattice.lattice_vec_a_comp(3,3)]);
        for i = 1:material_num
                plot_material(hax_comp, origin_point_set, org_idx{i}, Par_lattice.lattice_vec_a_comp, Par_lattice.lattice_constant_comp, Par_material.color_map{i}, 'computational');
        end
%         title('Material grids in computational domain')
    end        
end

function [ point_set_idx, reshaped_point_set ] = FLEEPS_Material_Locate( point_set, lattice_vec_a, Omega, Permutation, material_handle, periodic_judge )
    point_set_orig       = point_set*Omega;
    lattice_vec_a_orig_P = Omega'*lattice_vec_a;
    I = eye(3); P = I(:,Permutation); invP = P';
    invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
    lattice_vec_a_orig   = lattice_vec_a_orig_P(:,invPermutation);

    invAP  = inv(lattice_vec_a_orig_P);
    coef   = point_set_orig*invAP';

    shift_1 = floor(coef(:,1));
    shift_2 = floor(coef(:,2));
    shift_3 = floor(coef(:,3));

    point_set_orig = point_set_orig - [shift_1,shift_2,shift_3]*lattice_vec_a_orig_P'; 

    % Compute 27 shift values
    if strcmp(periodic_judge,'on') || strcmp(periodic_judge,'On')
        e = ones(3,1);  s = [-1;0;1];
        shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
        
        point_set_idx = cell(0);
        for i = 1:size(shift,1)
            shift_value    = shift(i,:)*lattice_vec_a_orig';
            point_set_idx_tmp  = material_handle(    point_set_orig(:,1)+shift_value(1),     point_set_orig(:,2)+shift_value(2),     point_set_orig(:,3)+shift_value(3), ...
                                                                lattice_vec_a_orig(:,1),                lattice_vec_a_orig(:,2),                lattice_vec_a_orig(:,3) );
            for j = 1:length(point_set_idx_tmp)
                if i == 1
                    point_set_idx{j} = [];
                end
                point_set_idx{j} = union( point_set_idx{j}, point_set_idx_tmp{j} );
            end
        end
    elseif strcmp(periodic_judge,'off') || strcmp(periodic_judge,'Off')
        point_set_idx  = material_handle(     point_set_orig(:,1),     point_set_orig(:,2),     point_set_orig(:,3), ...
                                          lattice_vec_a_orig(:,1), lattice_vec_a_orig(:,2), lattice_vec_a_orig(:,3) );
    end

    reshaped_point_set = point_set_orig;
end


function plotmesh3D(hax,mesh_x,mesh_y,mesh_z) 
    hold(hax, 'on')
    for i = 1:length(mesh_z)
        [tmp_grid_x,tmp_grid_y] = meshgrid(mesh_x,mesh_y);
        mesh(hax, tmp_grid_x, tmp_grid_y, mesh_z(i)*ones(size(tmp_grid_x)),'LineWidth',1 );
    end
    for i = 1:length(mesh_x)
        [tmp_grid_y,tmp_grid_z] = meshgrid(mesh_y,[mesh_z(1),mesh_z(end)] );
        mesh(hax, mesh_x(i)*ones(size(tmp_grid_y)), tmp_grid_y, tmp_grid_z,'LineWidth',1  );
    end
    hidden(hax, 'off');
    colormap(hax, [0,0,0]);
    axis(hax, 'equal');
%     hidden off 
%     colormap([0,0,0])
%     axis equal
%     axis off
end

function plot_material(hax,point_set, point_set_idx, lattice_vec_a, lattice_constant, color_map,lattice_vec_mode)
    plot3(hax,point_set(point_set_idx,1), point_set(point_set_idx,2), point_set(point_set_idx,3), 'o',...
          'MarkerSize',2,...
          'MarkeredgeColor',color_map,...
          'MarkerfaceColor',color_map);
    if isempty(lattice_vec_a) == 0
        FLEEPS_Plot_Parallelepiped(hax,[0,0,0],lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_constant,'color',lattice_vec_mode);
    end
    axis(hax, 'equal');
    axis(hax, 'tight');
%     axis off

    view(hax,[-15 0])
end