function [R1,M1] = filter_kmers_by_seq(S0,W0,min_overlap,max_ext)
% S0 = kmer sequences
% W0 = weight per kmer
%
% R = rows selected from original W matrix
% M = all rows merged with each selected row
% P = pwm for each selected row

if (nargin < 3)
    min_overlap = 4;
end
if (nargin < 4)
    max_ext = 3;
end

[~,i] = sortrows(W0);
S0 = S0(i);
W0 = W0(i,:);
n = size(S0,1);

% filter kmers with sequence overlaps
M1 = [];
m1 = [];
j = 1;
Sj = zeros(n,1);
for i = 1:size(S0,1)
    if (Sj(i) == 0)
        s1 = S0{i,1};
        k1 = i;
        
        % all kmers that are contained within s1
        k1 = [k1; i+find(cellfun(@isempty,regexp(s1,S0(i+1:end,1))) == 0)];

        % % all kmers with a sliding overlap with s1
        % N = kmer_overlap_all(S0(i+1:end,1),{s1},min_overlap,max_ext);
        % k1 = [k1; i+find(N==1)];

        % all kmers that contain s1
        if (length(s1) >= 4) 
            k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1)) == 0)];
        end

        % all kmers that contain s1 with 1 mismatch at the ends
        % if (length(s1) >= 5)
        %     k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1(1:end-1))) == 0)];
        %     k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1(2:end))) == 0)];
        % end

        % all kmers that contain s1 with 2 mismatches at the ends
        % if (length(s1) >= 6)
        %     k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1(1:end-2))) == 0)];
        %     k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1(3:end))) == 0)];
        %     k1 = [k1;i + find(cellfun(@isempty,regexp(S0(i+1:end,1),s1(2:end-1))) == 0)];
        % end

        k1 = setdiff(k1,find(Sj~=0));
        Sj(k1) = -1;
        Sj(i) = 1;
        M1{j} = [S0(k1,:) num2cell(W0(k1,:))];
        m1(j,1) = size(S0(k1,:),1);
        j = j+1;
    end
end
S1 = S0(Sj==1,1);
W1 = W0(Sj==1,:);
[R1,i] = sortrows([S1 num2cell([W1 m1])],2);
M1 = M1(i);
fprintf('K-mer filtering: %d input k-mers, %d final k-mers (%.0f%% filtered)\n',n,size(R1,1),100*(n-size(R1,1))/n);


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
