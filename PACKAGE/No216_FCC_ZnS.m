function material = No216_FCC_ZnS(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No216_FCC_ZnS
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cF8_216_c_a.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'ZnS';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cF8_216_c_a.html';

material.lattice_constant.a = 1;
[sphere_radius, cylinder_radius] = FLEEPS_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Mg
material.parameters{1}.name                 = 'Zn';
material.parameters{1}.color_map            = [66 99 145]/255;
material.parameters{1}.sphere_centers       = [    0 ,   0 ,    0 ];
                                                  
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [material.parameters{1}.sphere_centers;
                                                 0 , 0 , 1 ;
                                                 1 , 0 , 0 ;
                                                 0 , 1 , 0  ];
                                               
material.parameters{1}.cylinder_top_centers = 0.5* [  0.25 , 0.25 , 0.25 ;
                                                      0.25 , 0.25 , 1.25 ;
                                                      1.25 , 0.25 , 0.25 ;
                                                      0.25 , 1.25 , 0.25]; 
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 


%% Shape description for As
material.parameters{2}.name                 = 'S';
material.parameters{2}.color_map            = [255 242 0]/255;
material.parameters{2}.sphere_centers       = [    0.25 ,    0.25 ,     0.25 ];                                    
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ];
                                                                                                   
material.parameters{2}.cylinder_top_centers =  0.5*[  0.25 , 0.25 , 0.25 ;
                                                      0.25 , 0.25 , 1.25 ;
                                                      1.25 , 0.25 , 0.25 ;
                                                      0.25 , 1.25 , 0.25 ] ;
material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};