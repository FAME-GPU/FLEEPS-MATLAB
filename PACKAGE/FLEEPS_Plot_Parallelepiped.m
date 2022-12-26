function FLEEPS_Plot_Parallelepiped(hax,start_pt,vec1,vec2,vec3,lattice_constant,plot_mode, lattice_vec_mode)
    vec1     = reshape(vec1,3,1);
    vec2     = reshape(vec2,3,1);
    vec3     = reshape(vec3,3,1);
    start_pt = reshape(start_pt,3,1);
    if strcmp(plot_mode,'color') == 1
        position_a(:,1) = reshape((start_pt+vec1),3,1);
        position_a(:,2) = reshape((start_pt+vec2),3,1);
        position_a(:,3) = reshape((start_pt+vec3),3,1);
        
        rate = 1.1;
        if isfield(lattice_constant,'Permutation') == 1 && strcmp(lattice_vec_mode,'computational') == 1
%             plot3( [start_pt(1);position_a(1,1)], [start_pt(2);position_a(2,1)], [start_pt(3);position_a(3,1)], 'r-', 'LineWidth', 3);
%             plot3( [start_pt(1);position_a(1,2)], [start_pt(2);position_a(2,2)], [start_pt(3);position_a(3,2)], 'g-', 'LineWidth', 3);
%             plot3( [start_pt(1);position_a(1,3)], [start_pt(2);position_a(2,3)], [start_pt(3);position_a(3,3)], 'b-', 'LineWidth', 3);
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,1) position_a(2,1) position_a(3,1)], 'facealpha', 1, 'color', 'red', 'stemWidth', 0.02, 'tipWidth', 0.04); 
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,2) position_a(2,2) position_a(3,2)], 'facealpha', 1, 'color', 'green', 'stemWidth', 0.02, 'tipWidth', 0.04); 
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,3) position_a(2,3) position_a(3,3)], 'facealpha', 1, 'color', 'blue', 'stemWidth', 0.02, 'tipWidth', 0.04); 
            
            text(hax, rate*position_a(1,1), rate*position_a(2,1), rate*position_a(3,1), 'a_1', 'color', 'black', 'FontSize', 14);
            text(hax, rate*position_a(1,2), rate*position_a(2,2), rate*position_a(3,2), 'a_2', 'color', 'black', 'FontSize', 14);
            text(hax, rate*position_a(1,3), rate*position_a(2,3), rate*position_a(3,3), 'a_3', 'color', 'black', 'FontSize', 14);
        elseif strcmp(lattice_vec_mode,'original') == 1
%             position_a_perm = position_a(:,Permutation);
            
%             plot3( [start_pt(1);position_a_perm(1,1)], [start_pt(2);position_a_perm(2,1)], [start_pt(3);position_a_perm(3,1)], 'r-', 'LineWidth', 2);
%             plot3( [start_pt(1);position_a_perm(1,2)], [start_pt(2);position_a_perm(2,2)], [start_pt(3);position_a_perm(3,2)], 'g-', 'LineWidth', 2);
%             plot3( [start_pt(1);position_a_perm(1,3)], [start_pt(2);position_a_perm(2,3)], [start_pt(3);position_a_perm(3,3)], 'b-', 'LineWidth', 2);
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,1) position_a_perm(2,1) position_a_perm(3,1)], 'facealpha', 1, 'color', 'red', 'stemWidth', 0.02, 'tipWidth', 0.04); 
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,2) position_a_perm(2,2) position_a_perm(3,2)], 'facealpha', 1, 'color', 'green', 'stemWidth', 0.02, 'tipWidth', 0.04); 
%             mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,3) position_a_perm(2,3) position_a_perm(3,3)], 'facealpha', 1, 'color', 'blue', 'stemWidth', 0.02, 'tipWidth', 0.04); 
            
%             text(hax, rate*position_a(1,lattice_constant.Permutation(1)), rate*position_a(2,lattice_constant.Permutation(1)), rate*position_a(3,lattice_constant.Permutation(1)), '\bf{a}_1', 'color', 'black', 'FontSize', 14);
%             text(hax, rate*position_a(1,lattice_constant.Permutation(2)), rate*position_a(2,lattice_constant.Permutation(2)), rate*position_a(3,lattice_constant.Permutation(2)), '\bf{a}_2', 'color', 'black', 'FontSize', 14);
%             text(hax, rate*position_a(1,lattice_constant.Permutation(3)), rate*position_a(2,lattice_constant.Permutation(3)), rate*position_a(3,lattice_constant.Permutation(3)), '\bf{a}_3', 'color', 'black', 'FontSize', 14);
        end
    elseif strcmp(plot_mode,'no_color') == 1
        plot3(hax, [start_pt(1);start_pt(1)+vec1(1)'], [start_pt(2);start_pt(2)+vec1(2)'], [start_pt(3);start_pt(3)+vec1(3)'], 'k-',...
                   [start_pt(1);start_pt(1)+vec2(1)'], [start_pt(2);start_pt(2)+vec2(2)'], [start_pt(3);start_pt(3)+vec2(3)'], 'k-',...
                   [start_pt(1);start_pt(1)+vec3(1)'], [start_pt(2);start_pt(2)+vec3(2)'], [start_pt(3);start_pt(3)+vec3(3)'], 'k-');
    end
    vertex = [  start_pt(1) start_pt(2) start_pt(3);
                [start_pt(1) start_pt(2) start_pt(3)] + vec1';
                [start_pt(1) start_pt(2) start_pt(3)] + vec1'+vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1'+vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec2';];
    edge = [1 2;
            2 3;
            3 4;
            4 1;
            1 5;
            2 6;
            3 7;
            4 8;
            5 6;
            6 7;
            7 8;
            8 5];
    for i = 1:size(edge,1)        
        idx1 = edge(i,1);
        idx2 = edge(i,2);
        path = [vertex(idx1,:);vertex(idx2,:)];
        if i == 1
%             plot3( path(:,1), path(:,2), path(:,3), 'r-');
        elseif i == 4
%             plot3( path(:,1), path(:,2), path(:,3), 'g-');
        elseif i == 5
%             plot3( path(:,1), path(:,2), path(:,3), 'b-');
        else
            plot3(hax, path(:,1), path(:,2), path(:,3), 'k-')
        end
        hold(hax,'on')
%         if i == 1
%             plot3( path(1,1), path(1,2), path(1,3), 'ro')
%         end
    end

    end