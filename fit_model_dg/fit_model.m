function [P,E,X,R] = fit_model(t, x, mzt, k, LB, UB)
% t = times
% x = expression levels
% mzt = [max mzt time] [min offset time]
%
% P = <logX0,dg,t0,t1>

if ((nargin < 3) || isempty(mzt))
    mzt = [6 6];
end
if (nargin < 4)
    k = 10;
end
if (nargin < 6)
    LB = [];
    UB = [];
end

% dropouts
% calculate mean expression for each time + later samples
% use only samples with expression > (1/16)*mean
[~,~,w] = unique(t);
s = accumarray(w,x)./accumarray(w,1);
s = cumsum(s(end:-1:1))./(1:size(s,1))';
s = s(end:-1:1);
i = (x >= s(w)'-log2(16));
if (sum(i) < 4)
    [i;t;x;s(w)']'
end

% boundary
if (isempty(LB)||isempty(UB))
    mnT = max(min(t(i)),0);
    mxT = max(max(t(i)),6);

    %log(2)./[mzt(2) 0.03*mzt(2)]
    mxH = mzt(2);
    mnH = 0.03*mxH;

    LB = [-9 log(2)/mxH mnT mzt(2)];
    UB = [15 log(2)/mnH mzt(1) max(mxT+1,mzt(2))];
    if (UB(3) > mxT)
        LB(3) = mxT;
    end
    if (LB(4) > mxT)
        LB(4) = mxT+1;
        UB(4) = mxT+1;
    end
    
    nt = max(size(unique(t(i))));
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
X = nan(1,sum(i));
if (sum(i) < 4)
    fprintf('cannot perform fit (n=%d < 4)\n',sum(i));
    return;
end

% optimization
OPTS = optimset('TolX', 1e-5, 'TolFun', 1e-5, ...
    'MaxFunEvals',1e3, 'MaxIter', 1e3, 'display','off');

P0 = rnd_init(k, LB, UB);
for j = 1:k        
    %[p,e] = fminsearch(@errf, P0(j,:), OPTS, t(i), x(i));
    [p,e] = lsqnonlin(@errf, P0(j,:), LB, UB, OPTS, t(i), x(i));
    if (e < E)
        E = e;
        P = p;
    end
end
Xi = dg_eval_model(t(i),P);
R = corr(x(i)',Xi').^2;
X = dg_eval_model(t,P);

function e = errf(P, t, x)

w = dg_eval_model(t,P);
e = (w - x);
%e = sum(e.^2);

function P = rnd_init(k, LB, UB)

rng('shuffle');
n = size(LB,2);

P = zeros(k,n);
for i = 1:n
    P(:,i) = LB(i) + (UB(i) - LB(i))*rand(k,1);
end
