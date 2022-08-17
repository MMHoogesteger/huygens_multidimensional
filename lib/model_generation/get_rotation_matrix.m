function [A] = get_rotation_matrix(axis, angle)
%get_rotation_matrix returns rotation matrix for rotation about axis by
%angle
%   axis = 1,2 or 3, represents the axis of rotation;
%   angle, represents the angle of rotation
%   function returns matrix for a right-handed coordinate system
switch(axis)
    case 1
        A = [1 0 0 ; 0 cos(angle) sin(angle) ; 0 -sin(angle) cos(angle)];
    case 2
        A = [cos(angle) 0 -sin(angle) ; 0 1 0 ; sin(angle) 0 cos(angle)];
    case 3
        A = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1];
    otherwise
        error('no valid axis given')
end
end

