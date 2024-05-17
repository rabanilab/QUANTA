function p = standardize(data, r)
% data = each row is a gene
% r = # of significant digits

if (nargin < 2)
    r = NaN;
end
if (r > 20)
    r = 20;
end

% round expression values to the given significance level
if (~isnan(r))
    f = 10^r;
    d = round(data*f)/f;
else
    d = data;
end

% standardized data
[~,m] = size(d);
mn = repmat(mean_rows(d),1,m);
sd = repmat(std_rows(d),1,m);
p = (d-mn)./sd;



function mn = mean_rows(X)

[n,m] = size(X);
for(i=1:n)
    v = X(i,:);
    v = v(~isnan(v));
    if (isempty(v))
        mn(i,1) = 0;
    else
        mn(i,1) = mean(v);
    end
end
   

function sd = std_rows(X)

[n,m] = size(X);
for(i=1:n)
    v = X(i,:);
    v = v(~isnan(v));
    if (isempty(v))
        sd(i,1) = 1;
    else
        sd(i,1) = std(v);
    end
end
sd(sd == 0) = 1;
