function [R] = rot2(angle)
% rot2 returns rotation matrix for rotation of the system! by
% angle in x-y plane
% angle, represents the angle of rotation
% function returns matrix for a right-handed coordinate system, so positive
% angle direction is from x-axis to y-axis shortest angle.
R = [cos(angle) sin(angle); -sin(angle) cos(angle)];
end

