function [N,P,S] = pwm_score_all(PWM,seq)
% find and score best occurence of PWM in each sequence
%
% N = score; log2 P(seq|PWM)
% P = start position
% S = sequence

L = 'ACGT-';
for i = 1:4
    C(i) = sum([seq{:}]==L(i));
end
F = C./sum(C);

n = max(size(PWM));
m = size(seq,1);

N = zeros(m,n);
P = zeros(m,n);
S = cell(m,n);
for i = 1:n
    BW = repmat(F',1,size(PWM{i},2));
    [Ni,Pi,Si] = pwm_score(seq,PWM{i},BW,L);
    N(:,i) = Ni;
    P(:,i) = Pi;
    S(:,i) = Si;
end


function [opt_score,opt_pos,opt_seq] = pwm_score(seq,PWM,BW,L)
% calculate log-ratio score by given PWM matrix
% score = log2 P(seq|PWM)

nL = length(L);

% reassign PWMs with pseudo-counts
m = size(PWM,2);
log2pwm = add_pseudocounts(PWM);
log2bw = add_pseudocounts(BW);

% score sequences
k = size(seq,1);
opt_score = zeros(k,1);
opt_pos = zeros(k,1);
opt_seq = cell(k,1);
for i = 1:k
    seqi = [repmat('-',1,m) seq{i} repmat('-',1,m)];
    X = zeros(nL,length(seqi));
    for j = 1:nL
        X(j,seqi == L(j)) = 1;
    end
    PX = log2pwm'*X;
    sp = zeros(1,length(seqi)-m);
    for j = 1:m
        sp = sp + PX(j,j:end-m+(j-1));
    end
    BX = log2bw'*X;
    sb = zeros(1,length(seqi)-m);
    for j = 1:m
        sb = sb + BX(j,j:end-m+(j-1));
    end
    [max_score,max_pos] = max(sp-sb);
    opt_score(i) = max_score;
    opt_pos(i) = max_pos - m;
    opt_seq{i} = seqi(max_pos:(max_pos+m-1));
end

function log2pwm = add_pseudocounts(PWM)

pwmC = 20;
pseudoC = 1;

[n,m] = size(PWM);
C = round(pwmC.*PWM);
C = C + pseudoC;
C = [C; ones(1,m)];
log2pwm = log2(C) - repmat(log2(sum(C)),n+1,1);

