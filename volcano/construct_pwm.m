function [PWM,A,C] = construct_pwm(kmers,weights)
% generate a PWM matrix from input kmers by weights
%
% Input:
%  kmers = kmer sequences
%  weights = kmer weights
% Output:
%  N = PWM matrix
%  A = multiple sequence alignment
%  C = consensus sequence

if (nargin < 2)
    weights = ones(size(kmers));
end

L = {'A' 'C' 'G' 'T' '-'};
weightN = 1e3; % weighted PWM profile normalization factor
minI = 0.2; % minimal position information content to keep in final PWM
pseudoCNT = 25;

% multiple sequence alignment
if (size(kmers,1) == 2)
    A = multialign([kmers;kmers],'terminalGapAdjust',true,'GapOpen',100);
    A = A(1:2,:);
elseif (size(kmers,1) > 2)
    A = multialign(kmers,'terminalGapAdjust',true,'GapOpen',100);
else
    A = kmers{1};
end
m = size(A,2);

% weighted counts per position, PWM
PWM = zeros(5,m);
for j = 1:m
    [u,~,t] = unique(A(:,j));
    c = accumarray(t,weights);
    c = round(c*weightN);
    c(c<1) = 1;
    for k = 1:max(size(u))
        f = (strcmp(L,u(k)));
        PWM(f,j) = c(k);
    end
end

% add pseudo-counts, and calculate probabilities
PWM = PWM(1:4,:) + repmat(round(PWM(end,:)/4),4,1) + pseudoCNT;
PWM = PWM./repmat(sum(PWM,1),4,1);

% calculate information content
I = seqlogo(PWM,'DISPLAYLOGO','FALSE');
I = sum(I{2},1)';
k = find(I>=minI,1,'first'):find(I>=minI,1,'last');
PWM = PWM(:,k);
A = A(:,k);

% consensus
C = '';
for j = 1:size(PWM,2)
    [~,k] = max(PWM(:,j));
    C = [C L(k)];
end
C = [C{:}];
