function [ Point_idx ] = FLEEPS_Material_Locate_Handle_User_Defined( varargin )

X  = varargin{1}; Y  = varargin{2}; Z  = varargin{3};
a1 = varargin{4}; a2 = varargin{5}; a3 = varargin{6};
switch nargin
    case 9
        data_name       = varargin{7};
        sphere_radius   = varargin{8};
        cylinder_radius = varargin{9};
        Point_idx       = Locate_by_SphereCylinder(X,Y,Z,a1,a2,a3,data_name,sphere_radius,cylinder_radius);
    otherwise
        error('Too many or too few input parameters.');
end
end

function Point_idx = Locate_by_SphereCylinder(X,Y,Z,a1,a2,a3,data_name,sphere_radius,cylinder_radius)
Pmaterial = feval(data_name,sphere_radius,cylinder_radius);
Point_set = [X,Y,Z];
for i = 1:length(Pmaterial.parameters)
    Point_idx{i} = [];
    % find points inside spheres
    if isfield(Pmaterial.parameters{i},'sphere_centers') == 1
        for j = 1:size(Pmaterial.parameters{i}.sphere_centers,1)
            center = Pmaterial.parameters{i}.sphere_centers(j,:)*[a1,a2,a3]';
            Point_idx{i} = union(Point_idx{i},FLEEPS_Material_Locate_Sphere( Point_set, Pmaterial.parameters{i}.sphere_radius(j), center ));
        end
    end
    
    % find points inside cylinders
    if isfield(Pmaterial.parameters{i},'cylinder_bot_centers') == 1
        for j = 1:size(Pmaterial.parameters{i}.cylinder_bot_centers,1)
            bot_center = Pmaterial.parameters{i}.cylinder_bot_centers(j,:)*[a1,a2,a3]';
            top_center = Pmaterial.parameters{i}.cylinder_top_centers(j,:)*[a1,a2,a3]';
            Point_idx{i} =  union(Point_idx{i},FLEEPS_Material_Locate_Cylinder( Point_set, Pmaterial.parameters{i}.cylinder_radius(j), bot_center, top_center ));
        end
    end
end
end

function [ Point_idx ] = FLEEPS_Material_Locate_Sphere( Point_set, Sphere_radius, Sphere_center )
%% This function is used to determine whether a set of point in 3-dimensional euclidean space in seven spheres or not
% ==========================================================================================================
% Notice that the set of point need to be a m x 3 matrix, where m is the number of points
% The vector Sphere_center should be a 1 x 3 vector.
% The parameter for the location of sphere is setting in the function
%                        ' FLEEPS_Material_Parameter_Simple_Cube '
%
% ==========================================================================================================
% Sum = ( Point_set - Sphere_center ).^2 *ones(3,1);
Sum = ( Point_set(:,1) - Sphere_center(1) ).^2 + ( Point_set(:,2) - Sphere_center(2) ).^2 + ( Point_set(:,3) - Sphere_center(3) ).^2;
Point_idx =  find(Sum <= Sphere_radius*Sphere_radius);
end

function [ Point_idx ] = FLEEPS_Material_Locate_Cylinder( Point_set, Cylinder_radius, Cylinder_bot_center, Cylinder_top_center )
%% This function is used to determine whether a set of point in 3-dimensional euclidean space in a cylinder or not
% ==========================================================================================================
% Notice that the set of point need to be a m x 3 matrix, where m is the number of points
% The vectors Cylinder_bot_center and Cylinder_top_center should be two 1 x 3 vectors
% The parameter for the location of cylinder is setted in the function
%                        ' FLEEPS_Parameter_Simple_Cubic '
%
% ==========================================================================================================

if Cylinder_radius <= 0
    Point_idx = [];
else
    % Counting the number of points
    % Shifting the points
    Point_bot_set = Point_set - Cylinder_bot_center;
    Point_top_set = Point_set - Cylinder_top_center;
    % Computing the direction of cylinder
    Direct_vector = Cylinder_top_center - Cylinder_bot_center;
    % Normalizing the direction vector
    Direct_vector =       Direct_vector / norm( Direct_vector );
    % Constructing the array for vectorization computing
    %     Direct_vector_array = kron( ones( Point_num , 1 ) , Direct_vector );
    % Computing the distance for each point to the direction vector
    %     Point_line_dist   = cross( Point_bot_set , Direct_vector_array ).^2 * ones(3,1);
    Point_line_dist   =  ( Point_bot_set(:,2)*Direct_vector(3) - Point_bot_set(:,3)*Direct_vector(2) ).^2 + ...
        ( Point_bot_set(:,3)*Direct_vector(1) - Point_bot_set(:,1)*Direct_vector(3) ).^2 + ...
        ( Point_bot_set(:,1)*Direct_vector(2) - Point_bot_set(:,2)*Direct_vector(1) ).^2 ;
    % Computing norm of the projection vectors
    Point_line_bot_proj =    Point_bot_set(:,1)*Direct_vector(1) + Point_bot_set(:,2)*Direct_vector(2) + Point_bot_set(:,3)*Direct_vector(3);
    Point_line_top_proj = - (Point_top_set(:,1)*Direct_vector(1) + Point_top_set(:,2)*Direct_vector(2) + Point_top_set(:,3)*Direct_vector(3) );
    %     Point_line_proj     = min( [ Point_line_bot_proj , Point_line_top_proj ],[],2 );
    Point_line_proj     = min( Point_line_bot_proj,Point_line_top_proj );
    % Determining whether Point_set in a cylinder or not
    Point_idx_proj = find( Point_line_proj >= 0 );
    Point_idx_dist = find( Point_line_dist <= Cylinder_radius * Cylinder_radius );
    
    Point_idx = intersect(Point_idx_dist,Point_idx_proj);
end
end