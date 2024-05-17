function N = kldist_two(PWM1,PWM2,min_overlap)
% weight matrix comparison, using a method that has already been introduced by De Moor and colleagues (9,10)
% It is based on a symmetrized, position-averaged Kullback?Leibler distance.
% 
% To compare two weight matrices: 
% The shorter one is moved along the other. 
% All shifted positions that satisfy the following three conditions are considered:
%   1. >=50% of the shorter matrix overlaps with the longer matrix
%   2. overlap >= 4 positions 
%   3. the overlapping part of at least one matrix has to have a position-averaged entropy < 1 (natural logarithm). 
% For each such shift, a position-normalized dissimilarity score is calculated for the overlapping part.
% The smallest dissimilarity score is used to measure the overall similarity between the two matrices.
%
% Ref: PMC1160266

if (nargin < 3)
    min_overlap = 4;
end

n1 = size(PWM1,2);
n2 = size(PWM2,2);

if (n1>n2)
    N = kldist_two_ordered(PWM2,PWM1,min_overlap);
else
    N = kldist_two_ordered(PWM1,PWM2,min_overlap);
end

function N = kldist_two_ordered(P1,P2,min_overlap)

n1 = size(P1,2);
n2 = size(P2,2);
if (min(n1,n2) < min_overlap)
    min_overlap = min(n1,n2);
end

p = -1*max(floor(n1/2),n1-min_overlap) + 1;
N = inf;
while (p < n2)
    d = min(n1,n1 + p - 1);
    s1 = max(1-1*(p-1),1);
    e1 = min(s1+d,n1);
    s2 = max(p,1);
    e2 = min(s2+(e1-s1),n2);
    e1 = min(s1+(e2-s2),n1);
    l = e2-s2+1;
    %[[s1 e1 (e1-s1+1)]; ...
    %[s2 e2 (e2-s2+1)]]
    if (l >= min_overlap)
        P1i = P1(:,s1:e1)+eps;
        P2i = P2(:,s2:e2)+eps;
        ni = kldist(P1i,P2i);
        P1i = P1(:,[1:s1-1 e1+1:end])+eps;
        if (~isempty(P1i))
            ni = ni + kldist(P1i,0.25*ones(size(P1i)));
        end
        P2i = P2(:,[1:s2-1 e2+1:end])+eps;
        if (~isempty(P2i))
            ni = ni + kldist(P2i,0.25*ones(size(P2i)));
        end
        N = [N; ni];
    end
    p = p+1;
end
N = min(N);

function KL = kldist(P1,P2)

KL1 = sum(sum(P1.*log(P1./P2),2));
KL2 = sum(sum(P2.*log(P2./P1),2));
KL = (KL1+KL2)/2;

