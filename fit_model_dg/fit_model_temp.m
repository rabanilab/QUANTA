function [P,E,X,R] = fit_model_temp(t, x, timeB, k, LB, UB)
% t = times
% x = expression levels
% timeB = [max onset time] [min offset time]
%
% P = <logX0,dg,t0,t1>

if (nargin < 3)
    timeB = [];
end
if (nargin < 4)
    k = 10;
end
if (nargin < 6)
    LB = [];
    UB = [];
end

m = max(size(t));

% dropouts
% calculate mean expression for each time + later samples
% use only samples with expression > (1/16)*mean
i = [];
maxT = 0;
minT = inf;
for j = 1:m
    [~,~,w] = unique([t{j}]);
    s = accumarray(w,x{j})./accumarray(w,1);
    s = cumsum(s(end:-1:1))./(1:size(s,1))';
    s = s(end:-1:1);
    i{j} = (x{j} >= s(w)'-log2(16));
    maxT = max([t{j}(i{j}) maxT]);
    minT = min([t{j}(i{j}) minT]);
end

% timeB = [max onset time] [min offset time]
if (isempty(timeB))
    timeB(1) = 0.6*maxT;
    timeB(2) = 0.6*maxT;
end

% boundary
if (isempty(LB)||isempty(UB))
    mnT = max(minT,0);
    mxT = max(maxT,6);
    
    mxH = mxT;%timeB(2);
    mnH = 0.03*mxT;%0.03*mxH
    %DG: log(2)./[mxH mnH]

    LB = [-9 log(2)/mxH   mnT   timeB(2)];
    UB = [15 log(2)/mnH timeB(1) mxT+1];
    if (UB(3) > mxT)
        LB(3) = mxT;
    end
    if (LB(4) > mxT)
        LB(4) = mxT+1;
        UB(4) = mxT+1;
    end
    
    nt = max(size(unique([t{:}])));
    if (nt < 3)
        LB(3) = mnT;
        UB(3) = mnT;
    end
    if (nt < 4)
        LB(4) = mxT+1;
        UB(4) = mxT+1;
    end
end

% The Levenberg-Marquardt algorithm does not handle bound constraints
% and the trust-region-reflective algorithm requires at least as many
% equations as variables
P = nan(1,4);
E = inf;
R = nan;
X = nan(1,sum([i{:}]));
if (sum([i{:}]) < 4)
    fprintf('cannot perform fit (n=%d < 4)\n',sum([i{:}]));
    return;
end

% optimization
OPTS = optimset('TolX', 1e-5, 'TolFun', 1e-5, ...
    'MaxFunEvals',1e3, 'MaxIter', 1e3, 'display','off');

LB = [repmat(LB(1),1,m) LB(2:end)];
UB = [repmat(UB(1),1,m) UB(2:end)];
P0 = rnd_init(k, LB, UB);
for j = 1:k        
    %[p,e] = fminsearch(@errf, P0(j,:), OPTS, t(i), x(i));
    [p,e] = lsqnonlin(@errf, P0(j,:), LB, UB, OPTS, t,x,i);
    if (e < E)
        E = e;
        P = p;
    end
end

x1 = [];
x2 = [];
for j = 1:m
    y = dg_eval_model(t{j}(i{j}),[P(j) P(end-2:end)]);
    x1 = [x1 y];
    x2 = [x2 x{j}(i{j})];
end
R = corr(x1',x2').^2;
X = cell(1,m);
for j = 1:m
    X{j} = dg_eval_model(t{j},[P(j) P(end-2:end)]);
end

function e = errf(P, t,x,i)

m = max(size(t));
e = [];
for j = 1:m
    w = dg_eval_model(t{j}(i{j}),[P(j) P(end-2:end)]);
    e = [e (w - x{j}(i{j}))];
    %e = sum(e.^2);
end

function P = rnd_init(k, LB, UB)

rng('shuffle');
n = size(LB,2);

P = zeros(k,n);
for i = 1:n
    P(:,i) = LB(i) + (UB(i) - LB(i))*rand(k,1);
end
