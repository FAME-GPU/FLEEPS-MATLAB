function [ edge_center_x_set, edge_center_y_set, edge_center_z_set,...
           face_center_x_set, face_center_y_set, face_center_z_set,...
               origin_point_set ] = FLEEPS_Matrix_Grid( grid_nums, mesh_lens )
    n_grid_num = grid_nums(1)*grid_nums(2)*grid_nums(3);
    
    %% Create the Yee's scheme point
    % Creating the grids on the vertices
    x_grid  = (   0 : grid_nums(1)-  1 ).' *mesh_lens(1);
    y_grid  = (   0 : grid_nums(2)-  1 ).' *mesh_lens(2);
    z_grid  = (   0 : grid_nums(3)-  1 ).' *mesh_lens(3);
    % Creating the grids on the edge center
    x_shift_grid = ( 0.5 : grid_nums(1)-0.5 ).' *mesh_lens(1);
    y_shift_grid = ( 0.5 : grid_nums(2)-0.5 ).' *mesh_lens(2);
    z_shift_grid = ( 0.5 : grid_nums(3)-0.5 ).' *mesh_lens(3);
    %% Constructing points on edge for the electric field and on face for the magnetic field
    origin_point_set   = zeros(n_grid_num , 3);
    edge_center_x_set = zeros(n_grid_num , 3);
    edge_center_y_set = zeros(n_grid_num , 3);
    edge_center_z_set = zeros(n_grid_num , 3);
    face_center_x_set = zeros(n_grid_num , 3);
    face_center_y_set = zeros(n_grid_num , 3);
    face_center_z_set = zeros(n_grid_num , 3);
    % Constructing points on original grid
    origin_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );
    origin_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  );
    origin_point_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    % Constructing points on edge for the electric field
    edge_center_x_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) , x_shift_grid );
    edge_center_y_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );
    edge_center_z_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );

    edge_center_x_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  );
    edge_center_y_set(:,2) = kron( ones(grid_nums(3) , 1) , kron( y_shift_grid , ones(grid_nums(1) , 1)  )  );
    edge_center_z_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  ); 

    edge_center_x_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    edge_center_y_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    edge_center_z_set(:,3) = kron( z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    % Constructing points on face for the magnetic field
    face_center_x_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,      x_grid );
    face_center_y_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,x_shift_grid );
    face_center_z_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,x_shift_grid );

    face_center_x_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(y_shift_grid , ones(grid_nums(1) , 1)  )  );
    face_center_y_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(      y_grid , ones(grid_nums(1) , 1)  )  );
    face_center_z_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(y_shift_grid , ones(grid_nums(1) , 1)  )  ); 

    face_center_x_set(:,3) = kron(z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    face_center_y_set(:,3) = kron(z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    face_center_z_set(:,3) = kron(      z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
end