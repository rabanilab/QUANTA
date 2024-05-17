function [R1,P1,M1] = merge_kmers_by_seq(S,W,maxD,min_overlap,max_ext)
% S = kmer sequences
% W = weight per kmer
%
% R = rows selected from original W matrix
% M = all rows merged with each selected row
% P = pwm for each selected row

if (nargin < 3)
    maxD = 5;
end
if (nargin < 4)
    min_overlap = 4;
end
if (nargin < 5)
    max_ext = 3;
end

[~,i] = sortrows(W);
S = S(i);
W = W(i,:);
n = size(S,1);
fprintf('Input: %d k-mers\n',n);

% merge kmers with sequence overlaps
M1 = [];
P1 = [];
c1 = [];
m1 = [];
j = 1;
Sj = zeros(n,1);
for i = 1:size(S,1)
    if (Sj(i) == 0)
        s1 = S{i,1};
        k1 = i;
        
        % all kmers that are contained within s1
        k1 = [k1; i+find(cellfun(@isempty,regexp(s1,S(i+1:end,1))) == 0)];
        
        % all kmers with a sliding overlap with s1
        N = kmer_overlap_all(S(i+1:end,1),{s1},min_overlap,max_ext);
        k1 = [k1; i+find(N==1)];
        
        if (length(s1) >= 4) % all kmers that contain s1
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1)) == 0)];
        end
        if (length(s1) >= 5) % all kmers that contain s1 with 1 mismatch at the ends
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1(1:end-1))) == 0)];
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1(2:end))) == 0)];
        end
        if (length(s1) >= 6) % all kmers that contain s1 with 2 mismatches at the ends
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1(1:end-2))) == 0)];
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1(3:end))) == 0)];
            k1 = [k1;i + find(cellfun(@isempty,regexp(S(i+1:end,1),s1(2:end-1))) == 0)];
        end
        k1 = setdiff(k1,find(Sj~=0));
        Sj(k1) = -1;
        Sj(i) = 1;
        
        [PWM,~,CONS] = construct_pwm(S(k1,1),-1*log10(W(k1,1)));
        M1{j} = [S(k1,:) num2cell(W(k1,:))];
        m1(j,1) = size(S(k1,:),1);
        P1{j} = PWM;
        c1{j,1} = CONS;
        j = j+1;
    end
end
[R1,i] = sortrows([c1 S(Sj==1,1) num2cell([W(Sj==1,:) m1])],3);
M1 = M1(i);
P1 = P1(i);

n1 = size(R1,1);
n2 = inf;
i = 1;
while (n1 < n2)
    fprintf('iter %d: %d k-mers\n',i,n1);
    R2 = R1;
    M2 = M1;
    P2 = P1;
    n2 = n1;
    [R1,P1,M1] = merge_kmers_pwm(R2,M2,P2,n2,maxD,min_overlap);
    n1 = size(R1,1);
    i = i + 1;
end
fprintf('Final: %d k-mers\n',size(R1,1));


function [R1,P1,M1] = merge_kmers_pwm(R2,M2,P2,n2,maxD,min_overlap)

% calculate distance between matrices
N = pwm_compare(P2,min_overlap);

% merge pwms with minimal distance
k2 = [];
R1 = [];
P1 = [];
M1 = [];
for i = 1:n2
    if (ismember(i,k2))
        continue;
    end
    [m,mi] = min(N(i+1:end,i));
    mi = mi+i;
    
    % merge: select by information content
    if (m<maxD)
        k2 = [k2;mi];
        M = [M2{i}; M2{mi}];
        Mi = -1*log10(cell2mat(M(:,2)));
        [~,k] = sort(Mi,'descend');
        [p1,~,c1] = construct_pwm(M(k,1),Mi(k));
        R1 = [R1;[c1 R2(i,2:3) num2cell(size(Mi,1))]];
        M1 = [M1 {M(k,:)}];
        P1 = [P1 {p1}];
    else
        R1 = [R1;R2(i,:)];
        M1 = [M1 M2(i)];
        P1 = [P1 P2(i)];
    end
    k2 = [k2;i];
end

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


function N = kmer_overlap_all(k1,k2,min_overlap,max_ext)

m1 = size(k1,1);
m2 = size(k2,1);

N = zeros(m1,m2);
for i = 1:m1
    s1 = k1{i};
    n1 = length(s1);
    for j = 1:m2
        s2 = k2{j};
        n2 = length(s2);
        if (n1>n2)
            N(i,j) = kmer_overlap(s2,s1,min_overlap,max_ext);
        else
            N(i,j) = kmer_overlap(s1,s2,min_overlap,max_ext);
        end
    end
end

function N = kmer_overlap(k1,k2,min_overlap,max_ext)
% To compare two kmers: 
% The shorter one is moved along the other. 
% All shifted positions that satisfy the following two conditions are considered:
%   1. >=50% of the shorter kmer overlaps with the longer kmer
%   2. overlap >= 4 positions
%   3. non-overlap <= 3 positions
% For each such shift, we test the overall similarity between the two kmers

n1 = length(k1);
n2 = length(k2);

p = -1*max(floor(n1/2),n1-min_overlap) + 1;
N = 0;
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
    k1i = k1(:,s1:e1);
    k2i = k2(:,s2:e2);
    if ((l >= min_overlap) && (sum(k1i==k2i)>=l) && ((s1-1)+(s2-1)+(n1-e1)+(n2-e2)<=max_ext))
        N = 1;
        break;
    end
    p = p+1;
end
