function [P,E,X,R] = fit_model_2p(t, x, mzt, k, LB, UB)
% t = times
% x = expression levels
% mzt = [min mzt time] [max mzt time] [min offset time]
%
% P = <logX0,dg1,dg2,t0,t1>

if (nargin < 3)
    mzt = [3 6 6];
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
if (sum(i) < 5)
    [i;t;x;s(w)']'
    [LB;UB]
end

% boundary conditions
if (isempty(LB)||isempty(UB))
    mnT = max(min(t(i)),0);
    mnT = max(mnT,mzt(1));
    mxT = max(max(t(i)),6);

    %log(2)./[mzt(3) 0.03*mzt(3)]
    %log(2)./[-0.03*mzt(3) 0.03*mzt(3)]
    mxH = mzt(3);
    mnH = 0.03*mxH;

    LB = [-9 -1*log(2)/mnH log(2)/mxH  mnT  mzt(3)];
    UB = [15  log(2)/mnH log(2)/mnH mzt(2) mxT+1];
    if (UB(4) > mxT)
        UB(4) = mxT;
    end
    if (LB(5) > mxT)
        LB(5) = mxT+1;
        UB(5) = mxT+1;
    end
    
    nt = max(size(unique(t(i))));
    if (nt < 3)
        LB(4) = mnT;
        UB(4) = mnT;
        UB(2) = 0;
        LB(2) = 0;
    end
    if (nt < 5)
        LB(5) = mxT+1;
        UB(5) = mxT+1;
    end
end

% The Levenberg-Marquardt algorithm does not handle bound constraints
% and the trust-region-reflective algorithm requires at least as many
% equations as variables
P = nan(1,5);
E = inf;
R = nan;
X = nan(1,sum(i));
if (sum(i) < 5)
    fprintf('cannot perform fit (n=%d < 5)\n',sum(i));
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
Xi = dg_eval_model_2p(t(i),P);
R = corr(x(i)',Xi').^2;
X = dg_eval_model_2p(t,P);


function e = errf(P, t, x)

w = dg_eval_model_2p(t,P);
e = (w - x);
%e = sum(e.^2);

function P = rnd_init(k, LB, UB)

rng('shuffle');
n = size(LB,2);

P = zeros(k,n);
for i = 1:n
    P(:,i) = LB(i) + (UB(i) - LB(i))*rand(k,1);
end
