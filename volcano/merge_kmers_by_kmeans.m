function [PWM,Q,W,S] = merge_kmers_by_kmeans(pseq,pscore,iter,K,minK,maxK)
% k-means based approach to cluster kmers by sequence similarity
%  Centroids are defined by PWM
%  Distance is defined by PWM score
% 
% maxf = maximal fraction to include during initialization
% minK,maxK = k-mer lengths to test for initialization

if (nargin < 3)
    iter = 20;
end
if (nargin < 4)
    K = [];
end
if (nargin < 5)
    minK = 3;
end
if (nargin < 6)
    maxK = min(floor(size(pseq,1)/3),15);
end

if (isempty(K))
    K = init_select_k(pseq,pscore,minK,maxK);
    fprintf('*** Init K = %d [%d,%d] *** \n', K,minK,maxK);
else
    fprintf('*** Input K = %d *** \n', K);
end

S = 0;
PWM = [];
Q = [];
W = [];
if ((K>0)*(maxK>0) == 1)
    for i = 1:iter
        [PWMi,Qi,Wi,Si] = culster_kmers(pseq,pscore,K);
        fprintf('*** Iter %d (score = %.1f, %d motifs) *** \n', i, Si, size(PWMi,1));
        if ((Si > S) + (Si == S)*(size(PWMi,1)<size(PWM,1)) > 0)
            S = Si;
            PWM = PWMi;
            Q = Qi;
            W = Wi;
        end
    end
end
fprintf('*** Final (score = %.1f, %d motifs) *** \n', S, size(PWM,1));



function [PWM,Q,W,S] = culster_kmers(pseq,pscore,maxK)

% STEP 1: initialization
[PWM0,Q0] = init_pwm(pseq,pscore,maxK);
[xW0,xWid0,xPWM0,xQ0,xS0] = assign_to_pwm(pseq,pscore,pseq,PWM0,Q0);

% STEP 2: iterate
%   1. recalculate PWMs based on sequence assignments
%   2. assign all sequences to PWMs

% positive sequences
i1 = cell2mat(xQ0(:,3)) == 1;
j1 = cell2mat(xW0(:,4)) == 1;
[S1,PWM1,Q1,W1] = cluster_kmer_iterate(i1,j1,xW0,xWid0,xPWM0,xQ0,xS0);

% negative sequences
i2 = cell2mat(xQ0(:,3)) == -1;
j2 = cell2mat(xW0(:,4)) == -1;
[S2,PWM2,Q2,W2] = cluster_kmer_iterate(i2,j2,xW0,xWid0,xPWM0,xQ0,xS0);

% final
PWM = [PWM1;PWM2];
Q = [Q1;Q2];
W = [W1;W2];
S = sum(S1)+sum(S2);


function [S1,PWM1,Q1,W1] = cluster_kmer_iterate(i1,j1,xW0,xWid0,xPWM0,xQ0,xS0)

if (sum(j1) > 0)
    Wid1 = xWid0(j1);
    PWM1 = xPWM0(i1);
    Q1 = xQ0(i1,:);
    W1 = xW0(j1,:);
    S1 = xS0(i1,:);
    [W0,Wid0,PWM0,Q0,S0] = estimate_pwm(W1,Wid1,PWM1,Q1);
    
    % optimize
    k = 1;
    
    while (sum(S0)>sum(S1))
        S1 = S0;
        PWM1 = PWM0;
        Q1 = Q0;
        W1 = W0;
        Wid1 = Wid0;
        fprintf('* opt %d (%.1f, %d)\n',k,sum(S1),size(Q1,1));
        %[Q1 num2cell(S1)]
        
        % reassign to seeds
        [xW1,xWid1,xPWM1,xQ1] = assign_to_pwm(W1(:,1),cell2mat(W1(:,2)),W1(:,1),PWM1,Q1);
        [W0,Wid0,PWM0,Q0,S0] = estimate_pwm(xW1,xWid1,xPWM1,xQ1);
        k = k+1;
    end
    fprintf('* FINAL:\n');
    [Q1 num2cell(S1)]
else
    S1 = [];
    PWM1 = [];
    Q1 = [];
    W1 = [];
end


function [minK,maxK] = init_pwm_param()

minK = 4;
maxK = 6;

function K = init_select_k(pseq,pscore,minK,maxK)

iter = 10;
minC = 4;

% positive sequences
Kpos = 0;
k = pscore > 0;
if (sum(k)>0)
    Kpos = inf;
    for i = 1:iter
        Qpos = init_calc_pwm(pseq(k),pscore(k));
        W = cell2mat(Qpos(:,5));
        if (sum(W>=minC) < Kpos)
            Kpos = sum(W>=minC);
        end
    end
end

% negative sequences
Kneg = 0;
k = pscore < 0; 
if (sum(k)>0)
    Kneg = inf;
    for i = 1:iter
        Qneg = init_calc_pwm(pseq(k),pscore(k));
        W = cell2mat(Qneg(:,5));
        if (sum(W>=minC) < Kneg)
            Kneg = sum(W>=minC);
        end
    end
end

K = Kpos+Kneg;
if (K<minK)
    K = minK;
end
if (K>maxK)
    K = maxK;
end


function [PWM,Q] = init_pwm(pseq,pscore,maxK)

Q = [];

% positive sequences
k = pscore > 0;
if (sum(k)>0)
    [Qpos,PWMpos] = init_calc_pwm(pseq(k),pscore(k));
    Q = [Q; [Qpos(:,[2 5]) num2cell(ones(size(Qpos,1),1))]];
    fprintf('* INIT positive: %d sequences, %d selected seeds\n', sum(k), size(Qpos,1));
else
    PWMpos = [];
end

% negative sequences
k = pscore < 0; 
if (sum(k)>0)
    [Qneg,PWMneg] = init_calc_pwm(pseq(k),pscore(k));
    Q = [Q; [Qneg(:,[2 5]) num2cell(-1*ones(size(Qneg,1),1))]];
    fprintf('* INIT negative: %d sequences, %d selected seeds\n', sum(k), size(Qneg,1));
else
    PWMneg = [];
end
PWM = [PWMpos;PWMneg];

% take top PWMs
if(maxK<size(Q,1))
    Q = Q(1:maxK,:);
    PWM = PWM(1:maxK);
end

function [Q,PWM,L] = init_calc_pwm(pseq,pscore)

[minK,maxK] = init_pwm_param();

M = max(cellfun(@length,pseq));
if (M < minK)
    [Q,PWM,~,L] = init_select_seeds(pseq,pscore,M,M);
else
    [Q,PWM,~,L] = init_select_seeds(pseq,pscore,minK,maxK);
end


function [Q,PWM,A,seq_list] = init_select_seeds(pseq,pscore,mink,maxk)
% select seeds by score
% Q = [seed] [consensus] [seed length] [seed score] [number of peaks]

rng('shuffle');
[~,~,weightN,minI] = pwm_construct_param();

% count kmers in sequences (kmers x sequences)
Cid = [];
C = [];
plen = cellfun(@length,pseq);
for j = mink:maxk
    [kid,kc] = init_kmer_counts(pseq,j);
    Cid = [Cid;kid];
    C = [C;kc];
end

% weight[kmer,peak] = (fraction of kmer length out of overall peak length)*(peak score)
W = cellfun(@length,Cid)*abs(pscore./plen)';

% select seeds and optimize motifs
Kseed = [];
Kscore = [];
Kcnt = [];
Kcons = [];
PWM = [];
A = [];
seq_list = [];

% kmer score = sum of weights of all peaks that contain a kmer,
% normalized by kmer length relative to peak's length
I = (ones(1,size(C,2)) == 1);
kmerScore = sum((C(:,I)>0).*W(:,I),2);

while ((sum(I)>1)*(sum(kmerScore>0)>1) > 0)
        
    % select a random kmer
    %[s,j] = max(kmerScore);
    p = find(kmerScore>0);
    i = round(size(p,1)*rand(1));
    if (i <= 0)
        i = 1;
    end
    j = p(i);
    s = kmerScore(j);
    q = (C(j,:)>0).*(I>0) == 1;
    W(j,:) = 0;
    i = (sum((C>0)==repmat(q,size(C,1),1),2) == size(C,2));
    W(i,:) = 0;
    %fprintf('INIT seed = %s (score = %.2e) %d\n', Cid{j},s,sum(q));
    
    % calculate a consensus sequence
    [pwm,align,cons] = pwm_construct(pseq(q),abs(pscore(q)),Cid(j),[],weightN,minI);    
    
    Kseed = [Kseed;Cid(j)];
    Kscore = [Kscore;s];
    Kcons = [Kcons;{cons}];
    Kcnt = [Kcnt;sum(q>0)];
    PWM = [PWM;{pwm}];
    A = [A;{align}];
    seq_list = [seq_list;[num2cell(pscore(q)) pseq(q) repmat(Kcons(end),sum(q),1) repmat(Cid(j),sum(q),1)]];
        
    % keep only peaks that were not selected
    I = I.*(q==0) == 1;
    kmerScore = sum((C(:,I)>0).*W(:,I),2);
end
if (sum(I)>0)
    seq_list = [seq_list;[num2cell(pscore(I)) pseq(I) repmat({'NA'},sum(I),1) repmat({'NA'},sum(I),1)]];
end

Q = [Kseed Kcons num2cell([cellfun(@length,Kseed) Kscore Kcnt])];
[Q,i] = sortrows(Q,[-5 -4 -3]);
if (~isempty(PWM))
    PWM = PWM(i);
    A = A(i);
end


function [Cid,C] = init_kmer_counts(S,k)
% Count overlapping k-mers in input sequences

n = size(S,1);
Kseq = [];
Kid = [];
for j = 1:n
    s = S{j};
    l = length(s);
    for i = 1:l-k+1
        Kseq = [Kseq; s(i:i+k-1)];
        Kid = [Kid;j];
    end
end

if (~isempty(Kseq))
    Kseq = cellstr(Kseq);
    [Cid,~,t1] = unique(Kseq);
    W(:,1) = t1;
    W(:,2) = Kid;
    [W,~,t2] = unique(W,'rows');
    c = accumarray(t2,1);
    C = zeros(max(W(:,1)),n);
    for i = 1:size(c,1)
        C(W(i,1),W(i,2)) = c(i);
    end
else
    C = [];
    Cid = [];
end


function [pwmC,pseudoC,weightN,minI] = pwm_construct_param()

pwmC = 20;
pseudoC = 1;
weightN = 1e3; % weighted PWM profile normalization factor
minI = 0.1; % minimal position information content to keep in final PWM


function [N,A,C] = pwm_construct(kmers,weights,seed,A,weightN,minI)
% generate a PWM matrix from input kmers by weights
%
% Input:
%  kmers = kmer sequences
%  weights = kmer weights
%  seed = alignment seed (if exists)
%  weightN = weighted profile normalization factor
%  minI = minimal position information content to keep in final PWM
% Output:
%  N = PWM matrix
%  A = multiple sequence alignment
%  C = consensus sequence

% multiple sequence alignment
if (isempty(A))
    if (size(kmers,1)==2)
        A = multialign([kmers;kmers],'terminalGapAdjust',true,'GapOpen',100);
        A = A(1:2,:);
    elseif (~isempty(seed))
        A = pwm_seedalign(kmers,seed);
    elseif (size(kmers,1)>2)
        A = multialign(kmers,'terminalGapAdjust',true,'GapOpen',100);
    else
        A = kmers{1};
    end
end
[n,m] = size(A);

% weighted counts per position, PWM
L = {'A' 'C' 'G' 'T' '-'};
N = zeros(5,m);
for j = 1:m
    [u,~,t] = unique(A(:,j));
    c = accumarray(t,weights);
    c = round(c*weightN);
    c(c<1) = 1;
    for k = 1:max(size(u))
        f = (strcmp(L,u(k)));
        N(f,j) = c(k);
    end
end
N = N(1:4,:) + repmat(round(N(end,:)/4),4,1);
N = N./repmat(sum(N,1),4,1);

% calculate information content
I = seqlogo(N,'DISPLAYLOGO','FALSE');
I = sum(I{2},1)';
k = find(I>minI,1,'first'):find(I>minI,1,'last');
N = N(:,k);
A = A(:,k);

% consensus
C = '';
for j = 1:size(N,2)
    [~,k] = max(N(:,j));
    C = [C L(k)];
end
C = [C{:}];


function A = pwm_seedalign(kmers,seed)

l = cellfun(@length,kmers);
f = cellfun(@min,strfind(kmers,seed));
s = max(f)-f;
e = l+s;
e = max(e)-e;

n = size(kmers,1);
for i = 1:n
    A(i,:) = [repmat('-',1,s(i)) kmers{i,:} repmat('-',1,e(i))];
end


function [W,Wid,PWM,Q,S] = estimate_pwm(W0,Wid0,PWM0,Q0)

pseq = W0(:,1);
pscore = cell2mat(W0(:,2));
[~,~,weightN,minI] = pwm_construct_param();

Xid = unique(Wid0);
Xid = Xid(strcmp(Xid,'NA')==0);
PWMid = strcat(Q0(:,1),':',num2str(cell2mat(Q0(:,end))+1));
PWMid = regexprep(PWMid,':2',':P');
PWMid = regexprep(PWMid,':0',':N');

n = size(Xid,1);
PWM1 = cell(n,1);
Q1 = cell(n,3);
PWMw = zeros(n,1);
for i = 1:n
    j1 = strcmp(Wid0,Xid{i})==1;
    j2 = strcmp(PWMid,Xid{i})==1;
    if (sum(j2)>1)
        fprintf('%s [%d repeats]\n', Xid{i}, sum(j2));
        k = find(j2,1);
        j2 = zeros(size(j2));
        j2(k) = 1;
        j2 = (j2 == 1);
    end
    
    if (sum(j2)>0)
        [~,~,sseq] = pwm_score(W0(j1,1),PWM0{find(j2,1)});
        [PWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],cell2mat(sseq),weightN,minI);
    else
        [PWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],[],weightN,minI);
    end
    Q1(i,:) = [cons num2cell([sum(j1) 2*isempty(regexp(Xid{i},':N','once'))-1])];
    PWMw(i) = sum(sum(abs(cell2mat(W0(j1,2)))));
end

[W,Wid,PWM,Q,S] = assign_to_pwm(pseq,pscore,pseq,PWM1,Q1);
[S,i] = sort(S);
PWM = PWM(i);
Q = Q(i,:);


function [W,Wid,xPWM,xQ,xS] = assign_to_pwm(pseq,pscore,bgseq,PWM,Q)
% Output:
%  W = list of sequences [seq] [score] [motif score] [pos/neg] [assigned motif]
%  PWMid = motif ids
%  PWM = motif PWMs
%  PWMstat = peak statistics

minNA = 0;

% background distribution
[a,~,t] = unique([bgseq{:} 'ACGT']);
c = accumarray(t,1);
BG = c./sum(c);
%[cellstr(a') num2cell(BG)]

% assign peaks
PWMid = Q(:,1);
PWMs = cell2mat(Q(:,3));
%PWMn = cell2mat(PWMstat(:,2));

[W1,p1] = assign_seqs(BG,pseq(pscore>0),pscore(pscore>0),PWMs,PWMid,PWM,1);
%fprintf('Assign %d positive PWM models\n', sum(p1));
%[Q(p1,[1 3 2]) num2cell(100*PWMn(p1)./sum(PWMn))]

[W2,p2] = assign_seqs(BG,pseq(pscore<0),pscore(pscore<0),PWMs,PWMid,PWM,-1);
%fprintf('Assign %d negative PWM models\n', sum(p2));
%[Q(p2,[1 3 2]) num2cell(100*PWMn(p2)./sum(PWMn))]

W = [W1;W2];
W(cell2mat(W(:,3))<minNA,5) = {'NA'};
W(cell2mat(W(:,3))<minNA,3) = {0};
p = p1+p2 > 0;

Wid = strcat(W(:,end),':',num2str(cell2mat(W(:,end-1))+1));
Wid = regexprep(Wid,':2',':P');
Wid = regexprep(Wid,':0',':N');
Wid(strncmp(Wid,'NA:',3)) = {'NA'};
[u,~,t] = unique(Wid);
c = accumarray(t,1);
[u1,u2] = strtok(u,':');
xQ = [u1 num2cell([c 2*cellfun(@isempty,regexp(u2,':N'))-1])];
xS = accumarray(t,abs(cell2mat(W(:,3))));

xPWM = PWM(p);
xPWMid = strcat(PWMid(p),':',num2str(PWMs(p)+1));
xPWMid = regexprep(xPWMid,':2',':P');
xPWMid = regexprep(xPWMid,':0',':N');

[i,j] = ismember(xPWMid,u);
xPWM = xPWM(i);
xQ = xQ(j(j>0),:);
xS = xS(j(j>0));

function [W,p] = assign_seqs(BG,pseq,pscore,PWMs,PWMid,PWM,d)

n = size(pseq,1);
if (d>0)
    p = PWMs>0;
else
    p = PWMs<0;
end

W = [];
if (n>0)
    [S,N] = pwm_select(pseq,PWM(p),PWMid(p),BG);
    W = [W;[pseq num2cell([pscore N d*ones(n,1)]) S]];
end


function [ids,score,pos] = pwm_select(seq,PWM,pwm_ids,BG)
% assign sequences into best fitting PWM
% score = log2(P(seq|PWM)) - log2(P(seq|BG))

n = size(PWM,1);
S = ones(size(seq,1),n);
for i = 1:n
    l = size(PWM{i},2);
    bg = log2(pwm_score(seq,repmat(BG,1,l)));
    S(:,i) = log2(pwm_score(seq,PWM{i}));
    S(:,i) = S(:,i) - bg;
end

[score,pos] = max(S,[],2);
ids = pwm_ids(pos);


function [score,opt_pos,opt_seq] = pwm_score(seq,PWM)
% calculate log-ratio score by given PWM matrix
% score = P(seq|PWM)

% reassign PWMs with pseudo-counts
m = size(PWM,2);
log2pwm = add_pseudocounts(PWM);

% score sequences
k = size(seq,1);
score = zeros(k,1);
opt_pos = zeros(k,1);
opt_seq = cell(k,1);
for i = 1:k
    seqi = [repmat('-',1,m) seq{i} repmat('-',1,m)];
    max_score = -inf;
    max_pos = 0;
    for istart = 1:(length(seqi)-m)
        iend = istart + m - 1;
        scorei = get_score(seqi(istart:iend),log2pwm);
        if (scorei > max_score)
            max_score = scorei;
            max_pos = istart;
        end
    end
    score(i) = 2.^(max_score);
    opt_pos(i) = max_pos - m;
    opt_seq{i} = seqi(max_pos:(max_pos+m-1));
end

function score = get_score(seq,log2pwm)

L = 'ACGT-';
score = 0;
m = size(log2pwm,2);
for j = 1:m
    score = score + log2pwm(L==seq(j),j);
end

function log2pwm = add_pseudocounts(PWM)

[pwmC,pseudoC] = pwm_construct_param();

[n,m] = size(PWM);
C = round(pwmC.*PWM);
C = C + pseudoC;
C = [C; ones(1,m)];
log2pwm = log2(C) - repmat(log2(sum(C)),n+1,1);
