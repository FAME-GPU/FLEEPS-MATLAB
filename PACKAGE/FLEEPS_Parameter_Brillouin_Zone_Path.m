function Par_recip_lattice = FLEEPS_Parameter_Brillouin_Zone_Path( varargin )
    flag_usr = 0;
 
    switch nargin 
        case {0,1}
            error('To less input argument!')
        case 2
            wave_vec_array_usr = varargin{1};
            Par_recip_lattice  = varargin{2};
            flag_usr = 1;
            path_string = [];
            vertex = [];
            part_num = [];
        case 3
            part_num          = varargin{1};
            Par_lattice       = varargin{2};
            Par_recip_lattice = varargin{3};
            path_string       = default_path(Par_lattice.lattice_type, Par_lattice.lattice_constant_orig, Par_lattice.lattice_vec_a_orig);
        case 4
            part_num          = varargin{1};
            Par_lattice       = varargin{2};
            Par_recip_lattice = varargin{3};
            path_string       = varargin{4};
    end
    
    if flag_usr == 1
        wave_vec_array_usr = reshape(wave_vec_array_usr,3,size(wave_vec_array_usr,1)*size(wave_vec_array_usr,2)/3);
        wave_vec_array     = Par_recip_lattice.lattice_vec_b*wave_vec_array_usr;
    else
        [ vertex_orig ] = FLEEPS_Parameter_Brillouin_Zone_Point(Par_lattice.lattice_type, Par_lattice.lattice_constant_orig, Par_recip_lattice.lattice_vec_b_orig);
        vertex_str = fieldnames(vertex_orig);
        vertex = vertex_orig;
        for i = 1:length(vertex_str)
            eval( ['vertex.',vertex_str{i},' = Par_lattice.Omega*vertex_orig.',vertex_str{i},';'] );
        end
        
        n               = length(path_string);
        wave_vec_array  = [];
        if n == 1 
            eval( ['wave_vec_array = [ wave_vec_array, vertex.',path_string(1),'];']);
        end
        for i = 1:n-1
            if strcmp(path_string(i),'|') == 0 && strcmp(path_string(i+1),'|') == 0
                subpath_start_string = [ 'vertex.', path_string(i)   ];
                subpath_end_string   = [ 'vertex.', path_string(i+1) ];
                string_x = [ 'subpath(1,:) = linspace(',subpath_start_string,'(1),',subpath_end_string,'(1),', num2str(part_num),');'];
                string_y = [ 'subpath(2,:) = linspace(',subpath_start_string,'(2),',subpath_end_string,'(2),', num2str(part_num),');'];
                string_z = [ 'subpath(3,:) = linspace(',subpath_start_string,'(3),',subpath_end_string,'(3),', num2str(part_num),');'];

                eval(string_x);
                eval(string_y);
                eval(string_z);
                if i == n-1 || strcmp(path_string(i+2),'|') == 1
                    wave_vec_array = [wave_vec_array,subpath(:,1:end)];
                else
                    wave_vec_array = [wave_vec_array,subpath(:,1:end-1)];
                end
            end
        end
    end
    
    wave_vec_num           = size(wave_vec_array,2);

    Par_recip_lattice.wave_vec_array = wave_vec_array;
    Par_recip_lattice.wave_vec_num   = wave_vec_num; 
    Par_recip_lattice.path_string    = path_string;
    Par_recip_lattice.vertex         = vertex; 
    Par_recip_lattice.part_num       = part_num;
end

function path_string = default_path(lattice_type, lattice_constant, lattice_vec_a)
    switch lattice_type
        case 'simple_cubic'   
            path_string = 'GXMGRX'; %'GXMGRX|MR';
        case 'face_centered_cubic'
            path_string = 'GXWKGLUWLK'; %'GXWKGLUWLK|UX';
        case 'hexagonal'    
            path_string = 'GMKGALHA'; %'GMKGALHA|LM|KH';
    end
end

function [ vertex ] = FLEEPS_Parameter_Brillouin_Zone_Point(lattice_type, lattice_constant, reciprocal_lattice_vector_b)

switch lattice_type
    case 'simple_cubic'                                                      
        vertex.G = reciprocal_lattice_vector_b*[   0,   0,   0 ]';
        vertex.M = reciprocal_lattice_vector_b*[ 1/2, 1/2,   0 ]';
        vertex.R = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        vertex.X = reciprocal_lattice_vector_b*[   0, 1/2,   0 ]';

    case 'face_centered_cubic'                                                          
        vertex.G = reciprocal_lattice_vector_b*[ 0  , 0  , 0   ]';
        vertex.K = reciprocal_lattice_vector_b*[ 3/8, 3/8, 3/4 ]';
        vertex.L = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        vertex.U = reciprocal_lattice_vector_b*[ 5/8, 1/4, 5/8 ]';
        vertex.W = reciprocal_lattice_vector_b*[ 1/2, 1/4, 3/4 ]';
        vertex.X = reciprocal_lattice_vector_b*[ 1/2,   0, 1/2 ]';             

    case 'hexagonal'                                                                 
        vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
        vertex.A  = reciprocal_lattice_vector_b*[    0,   0, 1/2 ]';
        vertex.H  = reciprocal_lattice_vector_b*[  1/3, 1/3, 1/2 ]';
        vertex.K  = reciprocal_lattice_vector_b*[  1/3, 1/3,   0 ]';
        vertex.L  = reciprocal_lattice_vector_b*[  1/2,   0, 1/2 ]';
        vertex.M  = reciprocal_lattice_vector_b*[  1/2,   0,   0 ]';
end
end
