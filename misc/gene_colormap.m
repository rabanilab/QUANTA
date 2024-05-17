function cmp = gene_colormap(type)
% build a colormap for displaying genes:
% type = 
% 0: red,black,green
% 1: red,white,blue
% 2: orange,white,purple
% 3: orange,white,blue
% 4: red,white,orange
% 5: brown,white,green

if (nargin == 0)
    type = 0;
end

% parameters
dn0 = (1.0:-0.010:0.0)';
dn1 = (1.0:-0.005:0.5)';
dn2 = (1.0:-0.0075:0.25)';
up0 = (0.0: 0.010:1.0)';
up1 = (0.5: 0.005:1.0)';
up2 = (0.25:0.0075:1.0)';
s = round(size(dn0,1)*0.3);

% colormap
if (type == 0) % red,black,green
    low = zeros(size(up0));
    mid = zeros(s,1);
    
    red =   [low; mid; up0];
    green = [dn0; mid; low];
    blue =  [low; mid; low];
    
elseif (type == 1) % red,white,blue
    low = ones(size(up0));
    mid = ones(s,1);

    red =   [up0; mid; low];
    green = [up0; mid; dn0];
    blue =  [low; mid; dn0];
    
elseif (type == 2) % orange,white,purple
    low = ones(size(up0));
    mid = ones(s,1);
    
    red =   [up1; mid; low];
    green = [up0; mid; dn1];
    blue =  [low; mid; dn0];
    
elseif (type == 3) % orange,white,blue
    low = ones(size(up0));
    mid = ones(s,1);
    
    red =   [up0; mid; low];
    green = [up0; mid; dn1];
    blue =  [low; mid; dn0];

elseif (type == 4) % red,white,orange
    low = ones(size(up0));
    mid = ones(s,1);
    
    red =   [low; mid; low];
    green = [up1; mid; dn0];
    blue =  [up0; mid; dn0];

elseif (type == 5) % brown,white,green
    low = ones(size(up0));
    mid = ones(s,1);
    
    red =   [up0; mid; dn1];
    green = [up1; mid; dn2];
    blue =  [up0; mid; dn0];
end

cmp = [red green blue];
