function [m1,m2,m3,m4,t1,t2,t3,t4] = FLEEPS_Parameter_Boundary_Point(lattice_vec_a,lattice_constant,mesh_len, grid_num)
a1 = lattice_vec_a(:,1); a2 = lattice_vec_a(:,2); a3 = lattice_vec_a(:,3);

if lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2
    R3R8 = a2(1);
else
    R3R8 = a1(1) + a2(1);
end

% case 1(i-iv)
% if lattice_constant.theta_2 > 0  &&  lattice_constant.theta_2 < pi/2
if a3(1) > 0  &&  a3(2) > 0
    R9R1  = a3(1);
    R10R3 = a3(2);
    t1    = [0;0;0];
    t2    = -a1;
    %   case1-i
    if (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) <= a3(1)
        flag  = '1-i';
        R11R2 =  a3(1) - a2(1);
        t3    = -a1-a2; 
        t4    =  -a2;
    %   case1-ii
    elseif (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) >  a3(1)
        flag  = '1-ii';
        R11R2 = a1(1) + a3(1) - a2(1);
        t3    = -a2; 
        t4    =  a1-a2;
    %   case1-iii
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) <= a1(1)-a3(1)
        flag  = '1-iii';
        R11R2 = a3(1) - a2(1);
        t3    = -a1-a2; 
        t4    = -a2;
    %   case1-iv
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) >  a1(1)-a3(1)
        flag  = '1-iv';
        R11R2 = -a1(1) + a3(1) - a2(1);
        t3    = -2*a1-a2; 
        t4    = -a1-a2;
    end
end
% case 2(i-iv)
% if lattice_constant.theta_2 > pi/2  &&  lattice_constant.theta_2 < pi
if a3(1) < 0  &&  a3(2) > 0
    R9R1  = a1(1) + a3(1);
    R10R3 = a3(2);
    t1    = a1;
    t2    = [0;0;0];
    %   case2-i
    if (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) <= a1(1)+a3(1)
        flag  = '2-i';
        R11R2 = a1(1) - a2(1) + a3(1);
        t3    = -a2; 
        t4    =  a1-a2;
    %   case2-ii
    elseif (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) >  a1(1)+a3(1)
        flag  = '2-ii';
        R11R2 = 2*a1(1) - a2(1) + a3(1);
        t3    =  a1-a2; 
        t4    = 2*a1-a2;
    %   case2-iii
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) <= -a3(1)
        flag  = '2-iii';
        R11R2 = a1(1) + a3(1) - a2(1);
        t3    = -a2; 
        t4    =  a1-a2;
    %   case2-iv
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) >  -a3(1)
        flag  = '2-iv';
        R11R2 = - a2(1) + a3(1);
        t3    = -a1-a2; 
        t4    = -a2;
    end
end
% case 3(i-iv)
% if lattice_constant.theta_2 > pi  &&  lattice_constant.theta_2 < 3*pi/2
if a3(1) < 0  &&  a3(2) < 0
    R10R3 = a2(2) + a3(2);
    R11R2 = a1(1) + a3(1);
    t3    = [0;0;0];
    t4    = a1;
    %   case3-i
    if (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) <= -a3(1)
        flag  = '3-i';
        R9R1 = -a3(1) - a2(1);
        t1    = a1+a2;
        t2    = a2;
    %   case3-ii
    elseif (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) >  -a3(1)
        flag  = '3-ii';
        R9R1 =  a3(1) + a2(1);
        t1    =  a2;
        t2    = -a1+a2;
    %   case3-iii
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) <= a1(1)+a3(1)
        flag  = '3-iii';
        R9R1 =  a3(1) + a2(1) + a1(1);
        t1    = a1+a2;
        t2    = a2;
    %   case3-iv
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) >  a1(1)+a3(1)
        flag  = '3-iv';
        R9R1 =  a3(1) + a2(1) + 2*a1(1);
        t1    = 2*a1+a2;
        t2    = a1+a2;
    end
end
% case 4(i-iv)
% if lattice_constant.theta_2 > 3*pi/2  &&  lattice_constant.theta_2 < 2*pi
if a3(1) > 0  &&  a3(2) < 0
    R10R3 = a2(2) + a3(2);
    R11R2 = a3(1);
    t3    = -a1;
    t4    = [0;0;0];
    %   case4-i
    if (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) <= a1(1)-a3(1)
        flag  = '4-i';
        R9R1 =  a3(1) + a2(1);
        t1    = a2;
        t2    = -a1+a2;
    %   case4-ii
    elseif (lattice_constant.theta_3 > 0  &&  lattice_constant.theta_3 <= pi/2) && a2(1) >  a1(1)-a3(1)
        flag  = '4-ii';
        R9R1 =  a3(1) + a2(1) - a1(1);
        t1    = -a1+a2;
        t2    = -2*a1+a2;
    %   case4-iii
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) <= a3(1)
        flag  = '4-iii';
        R9R1 =  a3(1) + a2(1);
        t1    = a2;
        t2    = -a1+a2;
    %   case4-iv
    elseif (lattice_constant.theta_3 > pi/2  &&  lattice_constant.theta_3 < pi) && -a2(1) >  a3(1)
        flag  = '4-iv';
        R9R1 =  a3(1) + a2(1) + a1(1);
        t1    = a1+a2;
        t2    = a2;
    end
end
% flag
% R9R1
% R10R3
% R11R2
m2 = floor( (a3(1) + t1(1))/mesh_len(1)) ;
m3 = floor(R10R3/mesh_len(2)) ;


m1 = floor(R3R8/mesh_len(1));
% m2 = floor(R9R1/mesh_len(1)) ;
% m3 = floor(R10R3/mesh_len(2)) ;
% m4 = floor(R11R2/mesh_len(2))  ;
m4 = mod(m2 - m1, grid_num(1)) ;


if 0
    figure(1); cla
    hold on
    x = (0:grid_num(1) - 1)*mesh_len(1);
    y = (0:grid_num(2) - 1)*mesh_len(2);
    [X,Y] = meshgrid(x,y);
    X = X(:); Y = Y(:);
    plot(X,Y,'k.')

    plot([0,a1(1),a1(1),0,0],[0,0,a2(2),a2(2),0],'k-');
    plot([0,a1(1)],[0,0],'r-');
    plot([0,a2(1)],[0,a2(2)],'g-');
    if strcmp(flag,'2-i')
        plot([a1(1)-a2(1),a1(1)],[0,a2(2)],'g--');
    end
    if strcmp(flag,'4-iii')
        plot([a1(1),a1(1)+a2(1)],[0,a2(2)],'g--');
        plot([-a2(1),0],[0,a2(2)],'g--');
    end
    plot([0,a3(1)],[0,a3(2)],'b-');
    plot([a1(1)-R9R1,a1(1)-R9R1],[0,a2(2)-R10R3],'k-');
    plot([0,a1(1)],[a2(2)-R10R3,a2(2)-R10R3],'k-');
    plot([a1(1)-R11R2,a1(1)-R11R2],[a2(2)-R10R3,a2(2)],'k-');

    plot( (grid_num(1)-m1:grid_num(1)-1)*mesh_len(1),zeros(m1,1),'ro')
    plot( (grid_num(1)-m2:(grid_num(1)-1)) * mesh_len(1), zeros(m2,1),'g*')
    plot( zeros(m3,1), (grid_num(2)-m3:(grid_num(2)-1)) * mesh_len(2),'k*')
    plot( (grid_num(1)-m4:(grid_num(1)-1)) * mesh_len(1), a2(2)*ones(m4,1),'b*')
%     fprintf('m1 = %d, m2 = %d, m3 = %d, m4 = %d\n', m1,m2,m3,m4);
    axis equal
end