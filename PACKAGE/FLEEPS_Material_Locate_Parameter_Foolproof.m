function [sphere_radius, cylinder_radius] = FLEEPS_Material_Locate_Parameter_Foolproof(sphere_radius, cylinder_radius, material_num)
    if length(sphere_radius) < material_num
        sphere_radius = sphere_radius*ones(material_num,1);
    end
    if length(cylinder_radius) < material_num
        cylinder_radius = cylinder_radius*ones(material_num,1);
    end
end

function Pmaterial = FLEEPS_Material_Locate_Parameter_Get(file_name)
temp = 1;
eval(['Pmaterial=',file_name,'(temp,temp);']);
end



