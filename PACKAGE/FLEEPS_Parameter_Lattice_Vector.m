function [lattice_vec_a,lattice_vec_a_orig,length_a1,length_a2,length_a3,theta_1,theta_2,theta_3,Permutation] = FLEEPS_Parameter_Lattice_Vector(lattice_type,lattice_constant,lattice_vec_a)
% This routine return the lattice vectors and constants corresponding to
% given lattice type and lattice constants

%theta_1 is the angle between a1 and a2
%theta_2 is the angle between a1 and a3
%theta_3 is the angle between a2 and a3
    if exist('lattice_vec_a') ~= 1
        [a1_orig,a2_orig,a3_orig] = Lattice_vector_generate(lattice_type, lattice_constant);
    else
        a1_orig = lattice_vec_a(:,1);
        a2_orig = lattice_vec_a(:,2);
        a3_orig = lattice_vec_a(:,3);
    end
    
    
    lattice_vec_a_orig = [a1_orig,a2_orig,a3_orig];
    lattice_vec_a      = lattice_vec_a_orig;
    Permutation        = [1,2,3];
    
    length_a1_orig = norm(a1_orig);
    length_a2_orig = norm(a2_orig);
    length_a3_orig = norm(a3_orig);

    %% Condition 1: |a_1| > |a_2| & |a_3|
    [~,idx] = max([length_a1_orig,length_a2_orig,length_a3_orig]);
    lattice_vec_a(:,1)   = lattice_vec_a(:,idx);
    lattice_vec_a(:,idx) = lattice_vec_a_orig(:,1);
    Permutation([1,idx]) = Permutation([idx,1]);
    
    length_a1 = norm(lattice_vec_a(:,1));
    length_a2 = norm(lattice_vec_a(:,2));
    length_a3 = norm(lattice_vec_a(:,3));
    
    theta_1   = acos(lattice_vec_a(:,2)'*lattice_vec_a(:,3)/( length_a2*length_a3));
    theta_2   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,3)/( length_a1*length_a3));
    theta_3   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,2)/( length_a1*length_a2));
    
    %% Condition 2: ?a??sin£c_{£^}>((?a??)/(sin£c_{£^}))|cos£c_{£\}-cos£c_{£]}cos£c_{£^}|
    if length_a2*sin(theta_3) <= (length_a3/sin(theta_3))*abs(cos(theta_1)-cos(theta_2)*cos(theta_3))
        lattice_vec_a(:,[1,2,3]) = lattice_vec_a(:,[1,3,2]);
        Permutation(:,[1,2,3])   = Permutation(:,[1,3,2]);
        
        length_a1 = norm(lattice_vec_a(:,1));
        length_a2 = norm(lattice_vec_a(:,2));
        length_a3 = norm(lattice_vec_a(:,3));

        theta_1   = acos(lattice_vec_a(:,2)'*lattice_vec_a(:,3)/( length_a2*length_a3));
        theta_2   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,3)/( length_a1*length_a3));
        theta_3   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,2)/( length_a1*length_a2));
    end
   
    %% Final check
    if (length_a2 > length_a1) || (length_a3 > length_a1) || length_a2*sin(theta_3) <= (length_a3/sin(theta_3))*abs(cos(theta_1)-cos(theta_2)*cos(theta_3))
        error('The lattice constants does not suitable for computation!')
    end
    
    lattice_vec_a1 = [ length_a1; 0; 0 ];
    lattice_vec_a2 = [ length_a2*cos(theta_3); length_a2*sin(theta_3); 0];
    lattice_vec_a3 = [ length_a3*cos(theta_2); 
                       length_a3*( cos(theta_1) - cos(theta_3)*cos(theta_2) )/sin(theta_3);
                       length_a3*sqrt( 1 - cos(theta_3)^2 - cos(theta_2)^2 - cos(theta_1)^2 + ...
                                  2*cos(theta_3)*cos(theta_2)*cos(theta_1) )/sin(theta_3)];
    lattice_vec_a = [lattice_vec_a1, lattice_vec_a2, lattice_vec_a3];

end

function [a1,a2,a3] = Lattice_vector_generate(lattice_type, lattice_constant)
    %% Set lattice constants
    switch lattice_type
        %% Cubic system
        case 'simple_cubic'
            a    = lattice_constant.a;
            temp = a*eye(3);
            a1   = temp(:,1);
            a2   = temp(:,2);
            a3   = temp(:,3);

        case 'face_centered_cubic'
            a    = lattice_constant.a;
            temp = 0.5*a*[ 0, 1, 1  ;
                           1, 0, 1  ;
                           1, 1, 0 ];
            a1   = temp(:,1);
            a2   = temp(:,2);
            a3   = temp(:,3);     
        %% Hexagonal system
        case 'hexagonal' 
             a    = lattice_constant.a;
             c    = lattice_constant.c;
             %gamma = lattice_constant.gamma;
             temp = [ a,        -0.5*a, 0 ;
                      0, sqrt(3)*0.5*a, 0 ;
                      0,             0, c];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);

    end
end