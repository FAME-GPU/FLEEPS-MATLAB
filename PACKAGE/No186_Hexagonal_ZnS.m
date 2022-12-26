function material = No186_Hexagonal_ZnS(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No186_Hexagonal_ZnS
%
% This is the hexagonal analog of the zincblende lattice, i.e. the stacking of the ZnS 
% dimers is ABABAB¡K Replacing both the Zn and S atoms by C (or Si) gives the hexagonal 
% diamond structure. The ideal structure, where the nearest-neighbor environment of each
% atom is the same as in zincblende, is achieved when we take c/a=sqrt(8/3) and z2=1/8. 
% In the extreme case z2=1/2 this structure becomes the B_k (BN) structure. 
% Note that we have arbitrarily chosen the z1 parameter for the zinc atoms to be zero.
%        http://www.aflowlib.org/CrystalDatabase/AB_hP4_186_b_b.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
material.Prototype    = 'ZnS';
material.material_num = 2;

material.lattice_type = 'hexagonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_hP4_186_b_b.html';
material.z1           = 0.3748;
material.z2           = 0;

material.lattice_constant.a = 3.8227;
material.lattice_constant.c = sqrt(8/3)*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FLEEPS_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for S (Sphere)
material.parameters{1}.name                 = 'S';
material.parameters{1}.color_map            = [234 237 103]/255;

material.parameters{1}.sphere_centers       = [ 1/3, 2/3,     material.z1;
                                                2/3, 1/3, 0.5+material.z1 ];
                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Zn (Sphere)
material.parameters{2}.name                 = 'Zn';
material.parameters{2}.color_map            = [99 131 239]/255;

material.parameters{2}.sphere_centers       = [ 1/3, 2/3,     material.z2;
                                                2/3, 1/3, 0.5+material.z2];
                                           
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Shape description for S (Cylinder)
material.parameters{1}.cylinder_bot_centers =  [ material.parameters{1}.sphere_centers(1,:);
                                                 material.parameters{1}.sphere_centers(1,:);
                                                 material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(1,:);
                                                 material.parameters{1}.sphere_centers(1,:);
                                                 material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(2,:)];
material.parameters{1}.cylinder_top_centers =  .5*[ material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:);
                                                    material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:);
                                                    material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:);
                                                    material.parameters{2}.sphere_centers(1,:) + [0 0 1] + material.parameters{1}.sphere_centers(2,:)];
material.parameters{1}.cylinder_top_centers = [material.parameters{1}.cylinder_top_centers;
                                               material.parameters{1}.cylinder_top_centers(2,:) + 0.5*[0 1 0];
                                               material.parameters{1}.cylinder_top_centers(2,:) - 0.5*[1 0 0];
                                               material.parameters{1}.cylinder_top_centers(4,:) + 0.5*[1 0 0];
                                               material.parameters{1}.cylinder_top_centers(4,:) - 0.5*[0 1 0] ];                                                 
                                       
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 

%% Shape description for Zn (Cylinder)
material.parameters{2}.cylinder_bot_centers =  [ material.parameters{2}.sphere_centers(1,:);
                                                 material.parameters{2}.sphere_centers(2,:);
                                                 material.parameters{2}.sphere_centers(2,:);
                                                 material.parameters{2}.sphere_centers(1,:) + [0 0 1];
                                                 material.parameters{2}.sphere_centers(1,:) + [0 0 1];
                                                 material.parameters{2}.sphere_centers(1,:) + [0 0 1];
                                                 material.parameters{2}.sphere_centers(2,:);
                                                 material.parameters{2}.sphere_centers(2,:)];
material.parameters{2}.cylinder_top_centers =  .5*[ material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:);
                                                    material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:);
                                                    material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:);
                                                    material.parameters{2}.sphere_centers(1,:) + [0 0 1] + material.parameters{1}.sphere_centers(2,:)];
material.parameters{2}.cylinder_top_centers = [ material.parameters{2}.cylinder_top_centers;
                                                material.parameters{2}.cylinder_top_centers(4,:) + 0.5*[0 1 0];
                                                material.parameters{2}.cylinder_top_centers(4,:) - 0.5*[1 0 0];
                                                material.parameters{2}.cylinder_top_centers(2,:) + 0.5*[1 0 0];
                                                material.parameters{2}.cylinder_top_centers(2,:) - 0.5*[0 1 0];];                           
                                           
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};