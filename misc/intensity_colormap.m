function cmp = intensity_colormap(type)
% build a colormap for displaying intensities:
% type = {0: red; 1:green; 2:blue}

if (nargin == 0)
    type = 0;
end

% parameters
up  = [1:-0.001:0]';
s = size(up,1);

red = up;
green = up;
blue = up;

% colormap
if (type == 0)
    red = ones(s,1);  
elseif (type == 1)
    green = ones(s,1);
elseif (type == 2)
    blue = ones(s,1);
end

cmp = [red green blue];