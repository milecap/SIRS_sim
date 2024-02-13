function R = Rot(dir,a,d)

if strcmp(d,'lh')
    % clockwise rotations in a right handed system (my original version)
    % That means a left handed rotation on a right handed system
    if dir == 'x'
        R = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
    elseif dir == 'y'
        R = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
    elseif dir == 'z'
        R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
    else
        error('Error ocurred');
    end
elseif strcmp(d,'rh')
    % Counter clockwise rotations in a rigth handed system (reviewer version)
    % Right handed rotation on a right handed system
    if dir == 'x'
        R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    elseif dir == 'y'
        R = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
    elseif dir == 'z'
        R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    else
        error('Error ocurred');
    end
else
    error('Error ocurred');
end