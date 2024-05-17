function [d, i, c] = cluster_sort(D,k,n,minC)
% sort data by hierarchical clustering of the rows

if (nargin < 2)
    k = 1;
end
if (nargin < 3)
    n = 1;
end
if (nargin < 4)
    minC = 4;
end

% data standardization
R = standardize(D);

% sort matrix rows
m = size(R,1);
if (m <= 2)
    i = (1:m);
else
    pd = pdist(R, 'euclidean'); %, 'spearman'); % Spearman correlation
    l = linkage(pd, 'single'); % Shortest distance
    
    if (k == 2)
        a = clustergram(R, 'Cluster', 'column');
        close;
        i = str2double(get(a, 'RowLabels'))';
    else
        i = optimalleaforder(l, pd, 'Criteria', 'group', 'Transformation', 'linear')';%'adjacent'
    end
end
d = D(i,:);

% cut the tree
if ((n > 1) && (m > 2))
    c0 = cluster(l,'MaxClust',n);%'Cutoff',cutoff,'Criterion','distance'

    % re-number clusters
    x = unique(c0(i),'stable');
    y = unique(c0);
    c = zeros(size(c0));
    for j = 1:max(size(y))
        if (sum(c0==x(j)) >= minC)
            c(c0==x(j),1) = j;
        else
            c(c0==x(j),1) = 0;
        end
    end
    
    k = 0;
    for j = unique(c)'
        c(c==j) = -1*k;
        k = k+1;
    end
    c = -1*c;
    fprintf('%d clusters found\n', max(c));
else
    c = [];
end

