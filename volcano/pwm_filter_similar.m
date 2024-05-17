function [P2,S2] = pwm_filter_similar(P1,maxD,min_overlap)

if (nargin < 2)
    maxD = 5;
end
if (nargin < 3)
    min_overlap = 4;
end

P2 = [];
S2 = [];
n = max(size(P1));
fprintf('Input: %d PWMs\n',n);
if (n<2)
    return;
end

% iterative pairwise merge of PWMs
n1 = size(P1,2);
n2 = inf;
i = 1;
while (n2>n1)
    P2 = P1;
    n2 = n1;
    P1 = pwm_merge(P2,n2,maxD,min_overlap);
    n1 = size(P1,2);
    fprintf('iter %d: %d PWMs\n',i,n1);
    i = i+1;
end
fprintf('Final: %d PWMs\n',n1);

% calculate consensus sequences
for i = 1:n2
    S2{i,1} = pwm_cons(P2{i});
end


function P2 = pwm_merge(P1,n,maxD,min_overlap)

% calculate distance between matrices
N = pwm_compare(P1,min_overlap);

% merge pwms with minimal distance
k2 = [];
P2 = [];
for i = 1:n
    if (ismember(i,k2))
        continue;
    end
    [m,mi] = min(N(i+1:end,i));
    mi = mi+i;
    
    % merge: select by information content
    if (m<maxD)
        k2 = [k2;mi];
        i1 = seqlogo(P1{mi},'DISPLAYLOGO','FALSE');
        i1 = sum(i1{2},1);
        i2 = seqlogo(P1{i},'DISPLAYLOGO','FALSE');
        i2 = sum(i2{2},1);
        if (sum(i1) > sum(i2))
            P2 = [P2 P1(mi)];
        else
            P2 = [P2 P1(i)];
        end
    else
        P2 = [P2 P1(i)];
    end
    k2 = [k2;i];
end


function C = pwm_cons(PWM)

L = {'A' 'C' 'G' 'T'};
C = '';
for j = 1:size(PWM,2)
    [~,k] = max(PWM(:,j));
    C = [C L(k)];
end
C = [C{:}];

function N = pwm_compare(PWM,min_overlap)
% compare a list of PWM matrices with itself

n = max(size(PWM));
N = zeros(n,n);

for i = 1:n
    for j = i+1:n
        N(i,j) = kldist_two(PWM{i},PWM{j},min_overlap);
        N(j,i) = N(i,j);
    end
end

function N = pwm_compare_two(PWM1,PWM2)
% compare two lists of PWM matrices

m1 = max(size(PWM1));
m2 = max(size(PWM2));
N = zeros(m1,m2);

for i = 1:m1
    for j = 1:m2
        N(i,j) = kldist_two(PWM1{i},PWM2{j});
    end
end
