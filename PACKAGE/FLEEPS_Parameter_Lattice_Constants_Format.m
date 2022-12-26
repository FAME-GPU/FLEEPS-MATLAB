function lattice_constant = FLEEPS_Parameter_Lattice_Constants_Format( varargin )
    lattice_type = varargin{1};
    if isstruct(varargin{2}) == 1
        lattice_constant = varargin{2};
        lattice_constant = check_lattice_constant( lattice_type, lattice_constant );
    else
        lattice_vec_a = varargin{2};
        % lattice_constant = check_lattice_vec( lattice_type, lattice_vec_a );
    end

end

function lattice_constant = check_lattice_constant( lattice_type, lattice_constant )

% This routine reture pretreatment lattice constants
    lattice_constant.old = lattice_constant;
    switch lattice_type  
        case {'simple_cubic','face_centered_cubic'}
              %% Cubic system
        % In cubic system, lattice constant must satisfying
        %               a = b = c, £\=£]=£^=£k/2
            if isfield(lattice_constant,{'a','b','c'}) ~= 1
                error('There has only one lattice constant ''a'' should be specified. Notice that, in cubic system, lattice constant must satisfying a = b = c, £\=£]=£^=£k/2');
            end

            lattice_constant.b     = lattice_constant.a;
            lattice_constant.c     = lattice_constant.a;
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'hexagonal'}
              %% Hexagonal systems
        % In hexagonal and rhombohedral systems, lattice constant must satisfying
        %               a = b, £\=£]=£k/2, £^=2£k/3
            if sum(isfield(lattice_constant,{'a','c'})) ~= 2
                error('The lattice constant ''a'' or ''c'' were not specified. Notice that, in hexagonal and rhombohedral systems, lattice constant must satisfying a = b, £\=£]=£k/2, £^=2£k/3');
            end
            lattice_constant.b     = lattice_constant.a;
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = 2*pi/3;
    end
end