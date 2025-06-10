function [maxQ,maxPWM] = fit_pwms(PEAKs,max_peaks)
% Input:
%  PEAKs = list of peaks [weight][sequence]
% Output:
%  Q = [consensus] [number of peaks] [POS/NEG]
%  PWM = position weight matrix per peak motif

if (nargin < 2)
    max_peaks = 2000;
end

% seeds
min_seed_kmer = 4;
max_seed_kmer = 8;
max_seeds = 15;

% pwm fit
minK = 1;
maxK = 8;
maxK_iter = 10;%25;
max_iter = 25;
score_type = 2; % 1=enhanced BIC, 2=BIC, 3=AIC

% minI = minimal position information content to keep in final PWM
minI = 0.5;


% ------------------------
% limit the number of peaks, if too high
% ------------------------
n = size(PEAKs,1);
if (n > max_peaks)
    PEAKs = PEAKs(1:max_peaks,:);
end
fprintf('[1] Peaks: n=%d (%.0f%%)\n',size(PEAKs,1),100*size(PEAKs,1)/n);

% ------------------------
% initialize PWMs
% ------------------------
bgseq = PEAKs(:,2);
pseq = PEAKs(:,2);
pscore = cell2mat(PEAKs(:,1));
[Q0,PWM0] = select_seeds(pseq,pscore,min_seed_kmer,max_seed_kmer);

Q0 = [Q0(:,[2 5]) num2cell(ones(size(Q0,1),1))];
n = size(Q0,1);
j = min(size(Q0,1),max_seeds);
Q0 = Q0(1:j,:);
PWM0 = PWM0(1:j);

% remove from PWMs positions with low information content
for i = 1:max(size(PWM0))
    N = PWM0{i};
    A = Q0{i,1};
    I = seqlogo(N,'DISPLAYLOGO','FALSE');
    I = sum(I{2},1)';
    k = find(I>=minI,1,'first'):find(I>=minI,1,'last');
    if (isempty(k))
        k = find(I>=0.5*max(I),1,'first'):find(I>=0.5*max(I),1,'last');
    end
    PWM0{i} = N(:,k);
    Q0{i,1} = A(k);
end

fprintf('[1] Seeds: %d seeds\n', n);
fprintf('[2] Seeds: top %d seeds selected (%.1f%%):\n', j, 100*j/n);
Q0

% ------------------------
% fit PWMs to peaks
% ------------------------
% Iteratively:
%  [1] calculate a new PWM matrix for all assigned peaks
%  [2] assign peaks to new PWMs
if (maxK > size(Q0,1))
    maxK = size(Q0,1);
end
fprintf('K range: [%d,%d] PWMs to fit\n',minK,maxK);

maxS = zeros(1,maxK);
sPWM = cell(1,maxK);
sQ = cell(1,maxK);
parfor ik = minK:maxK
    [kPWM,kQ,~,kS] = peaks_fit_pwm_k(PWM0,Q0,pseq,pscore,bgseq,ik,score_type,maxK_iter,max_iter);
    maxS(ik) = kS;
    sPWM{ik} = kPWM;
    sQ{ik} = kQ;
end

num2cell([(1:maxK)' maxS'])
[~,i] = max(maxS);
maxPWM = sPWM{i};
maxQ = sQ{i};

% ------------------------
% remove from PWMs positions with low information content
% ------------------------
rm = 1;
rmi = 0;
while ((rm>0) && (rmi<max_iter))
    rm = 0;
    for i = 1:max(size(maxPWM))
        N = maxPWM{i};
        A = maxQ{i,1};
        I = seqlogo(N,'DISPLAYLOGO','FALSE');
        I = sum(I{2},1)';
        k = find(I>=minI,1,'first'):find(I>=minI,1,'last');
        if (isempty(k))
            k = find(I>=0.5*max(I),1,'first'):find(I>=0.5*max(I),1,'last');
            rm = 1;
            rmi = rmi+1;
        end
        maxPWM{i} = N(:,k);
        maxQ{i,1} = A(k);
    end
    [W0,Wid0,xPWM0,xPWMstat0] = peaks_assign(pseq,pscore,pseq,maxPWM,maxQ(:,1));
    [maxPWM,maxQ] = peaks_fit_pwm(W0,Wid0,xPWM0,xPWMstat0);
end

[~,i] = sortrows(maxQ(:,2),-1);
maxQ = maxQ(i,1:2);
maxPWM = maxPWM(i);
q = cell2mat(maxQ(:,2));
maxQ(:,3) = num2cell(q./sum(q));

fprintf('FINAL assigned peaks:\n');
maxQ


% --------------------------------------------------------------------
% select seeds
% --------------------------------------------------------------------

function [Q,PWM,A,seq_list] = select_seeds(pseq,pscore,mink,maxk)
% select seeds by score
% Q = [seed] [consensus] [seed length] [seed score] [number of peaks]

% count kmers in peak sequences
Cid = [];
C = [];
plen = cellfun(@length,pseq);
for j = mink:maxk
    [kid,kc] = kmer_counts(pseq,j);
    Cid = [Cid;kid];
    C = [C;kc];
end

% weight[seed,peak] = (fraction of seed length out of overall peak length)*(peak score)
W = cellfun(@length,Cid)*abs(pscore./plen)';

% select seeds and optimize motifs
Kseed = [];
Kscore = [];
Kcnt = [];
Kcons = [];
PWM = [];
A = [];
seq_list = [];

% seed score = sum of weights of all peaks that contain a seed,
%              normalized by seed length relative to peak's length
I = (ones(1,size(C,2)) == 1);
seedScore = sum((C(:,I)>0).*W(:,I),2);

while ((sum(I)>0)*(sum(seedScore>0)>0) > 0)
        
    % select a kmer with maximal score
    [s,j] = max(seedScore);
    q = (C(j,:)>0).*(I>0) == 1;
    W(j,:) = 0;
    i = (sum((C>0)==repmat(q,size(C,1),1),2) == size(C,2));
    W(i,:) = 0;
    fprintf('select seed: %s (score = %.2e) %d\n', Cid{j},s,sum(q));
    
    % calculate a consensus sequence
    [pwm,align,cons] = pwm_construct(pseq(q),abs(pscore(q)),Cid(j),[]);        
    Kseed = [Kseed;Cid(j)];
    Kscore = [Kscore;s];
    Kcons = [Kcons;{cons}];
    Kcnt = [Kcnt;sum(q>0)];
    PWM = [PWM;{pwm}];
    A = [A;{align}];
    seq_list = [seq_list;[num2cell(pscore(q)) pseq(q) repmat(Kcons(end),sum(q),1) repmat(Cid(j),sum(q),1)]];
        
    % keep only seeds that were not already selected
    I = I.*(q==0) == 1;
    seedScore = sum((C(:,I)>0).*W(:,I),2);
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

function [Cid,C] = kmer_counts(S,k)
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


% --------------------------------------------------------------------
% fit PWMs
% --------------------------------------------------------------------
function [kPWM,kQ,kS,kSL] = peaks_fit_pwm_k(PWM0,Q0,pseq,pscore,bgseq,ik,score_type,maxK_iter,max_iter)

kSL = -inf;
kS = -inf;
kPWM = [];
kQ = [];
iiter = min(ceil(0.75*nchoosek(size(PWM0,1),ik)),maxK_iter);
fprintf('K = %d PWMs to fit (%d iters)\n',ik,iiter);

for it = 1:iiter
    I = randperm(size(PWM0,1))';
    I(I>ik) = 0;
    [W0,Wid0,xPWM0,xPWMstat0,xPWMw0] = peaks_assign(pseq,pscore,bgseq,PWM0(I>0),Q0(I>0,1));

    % Iteratively:
    %  [1] calculate a new PWM matrix for all assigned peaks
    %  [2] assign peaks to new PWMs
    Si = -inf;
    S1 = sum(cell2mat(W0(:,3)));
    iter = 0;
    while ((S1>Si)*(iter<max_iter) == 1)
        PWM = xPWM0;
        Q = xPWMstat0;
        pwm_w = xPWMw0;
        Si = S1;

        % calculate new PWM matrices for all assigned peaks
        %fprintf('[iter %d] %.2e (n=%d)\n', iter, S1, size(PWM0,1));
        [PWM1,PWMstat1] = peaks_fit_pwm(W0,Wid0,xPWM0,xPWMstat0);

        % assign peaks to PWMs
        [W0,Wid0,xPWM0,xPWMstat0,xPWMw0] = peaks_assign(pseq,pscore,pseq,PWM1,PWMstat1(:,1));
        %fprintf('[iter %d] assigned peaks:\n', iter);
        %[PWMstat0 num2cell(PWMw0)]

        S1 = sum(cell2mat(W0(:,3)));
        iter = iter + 1;
    end

    % sort motifs by score
    [~,i] = sort(pwm_w);
    PWM = PWM(i);
    Q = Q(i,:);

    % scoring the result
    N = (size([PWM{:}],1)-1)*size([PWM{:}],2);
    BIC1 = 2*log(2)*Si - 10*log(length(pseq))*N; % enhanced BIC score 
    BIC2 = 2*log(2)*Si - log(length([pseq{:}]))*N; % BIC score
    AIC = 2*(log(2)*Si - N); % AIC score
    fprintf('final K=%d i=%d logL=%.1f (%d param, eBIC %.0f, BIC %.0f, AIC %.0f)\n\n', ...
        ik,it,Si,N,BIC1,BIC2,AIC);
    if (Si > kS)
        kS = Si;
        kPWM = PWM;
        kQ = Q;
        if (score_type == 3)
            kSL = AIC;
        elseif (score_type == 2)
            kSL = BIC2;
        else
            kSL = BIC1;
        end
    end
end

function [W,Wid,xPWM,xPWMstat,xPWMw] = peaks_assign(pseq,pscore,bgseq,PWM,PWMid)
% assign peaks into best fitting PWM
%   score = log2(P(seq|PWM)) - log2(P(seq|BG))
%
% Output:
%  W = list of peaks [peak seq] [peak score] [motif score] [pos/neg] [assigned motif]
%  PWMid = motif ids
%  PWM = motif PWMs
%  PWMstat = peak statistics [peak seq] [peak count] [pos/neg]

% select a PWM for each peak
n = size(pseq,1);
W = [];
if (n>0)
    S = pwm_score_all(pseq,bgseq,PWM);
    [N,pos] = max(S,[],2);
    S = PWMid(pos);
    %S(N<0) = {'NA'};
    W = [W;[pseq num2cell([pscore N ones(n,1)]) S]];
end

Wid = W(:,end);
[u,~,t] = unique(Wid);
c = accumarray(t,1);
w = accumarray(t,abs(cell2mat(W(:,2))));
xPWMstat = [u num2cell([c ones(size(w))])];
xPWMw = w;
xPWM = PWM;
xPWMid = PWMid;
[i,j] = ismember(xPWMid,u);
xPWM = xPWM(i);
xPWMstat = xPWMstat(j(j>0),:);
xPWMw = xPWMw(j(j>0));


function [PWM1,PWMstat1] = peaks_fit_pwm(W0,Wid0,PWM0,PWMstat0)
% calculate a new PWM matrix for all assigned peaks

Xid = unique(Wid0);
n = size(Xid,1);

id = PWMstat0(:,1);
PWM1 = cell(n,1);
PWMstat1 = cell(n,3);
W1 = zeros(n,1);
for i = 1:n
    j1 = strcmp(Wid0,Xid{i})==1;
    j2 = strcmp(id,Xid{i})==1;

    % found in the list of existing PWMs
    if (sum(j2)>0) 
        k = find(j2,1,'first');
        [~,~,sseq] = pwm_score(W0(j1,1),PWM0{k});
        [PWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],cell2mat(sseq));
    % % generate a new PWM
    % else 
    %     [PWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],[]);
    end
    PWMstat1(i,:) = [cons num2cell([sum(j1) 1])];
    W1(i) = sum(sum(abs(cell2mat(W0(j1,2)))));
    %fprintf('[%d] %s -> %s (w = %.2f)\n', i,Xid{i},cons,W1(i));
end

function S = pwm_score_all(pseq,bgseq,PWM)
% S = score matrix of size [pseq] x [PWMs]
% score = log2 P(seq|PWM)/P(seq/BG)

% calculate a background distribution
bg_pseudo = repmat('ACGT',1,10);
[a,~,t] = unique([bgseq{:} bg_pseudo]);
c = accumarray(t,1);
BG = c./sum(c);
%[cellstr(a') num2cell(BG)]

% calculate scores
m = max(size(PWM));
S = ones(size(pseq,1),m);
for i = 1:m
    l = size(PWM{i},2);
    bg = log2(pwm_score(pseq,repmat(BG,1,l)));
    S(:,i) = log2(pwm_score(pseq,PWM{i}));
    S(:,i) = S(:,i) - bg;
end

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

function [pwmC,pseudoC,weightN] = pwm_param()

pwmC = 20;
pseudoC = 1;
weightN = 1e3;

function log2pwm = add_pseudocounts(PWM)

% pwmC = pwm counts
% pseudoC = pseudo-counts
[pwmC,pseudoC] = pwm_param();

[n,m] = size(PWM);
C = round(pwmC.*PWM);
C = C + pseudoC;
C = [C; pseudoC*ones(1,m)];
log2pwm = log2(C) - repmat(log2(sum(C)),n+1,1);

function [N,A,C] = pwm_construct(kmers,weights,seed,A)
% generate a PWM matrix from input kmers by weights
%
% Input:
%  kmers = kmer sequences
%  weights = kmer weights
%  seed = alignment seed (if exists)
%  A = alignment
% Output:
%  N = PWM matrix
%  A = multiple sequence alignment
%  C = consensus sequence

%  weightN = weighted profile normalization factor
[~,~,weightN] = pwm_param();

% multiple sequence alignment
if (isempty(A))
    if (size(kmers,1)==1)
        A = kmers{1};
    elseif (size(kmers,1)==2)
        A = multialign([kmers;kmers],'terminalGapAdjust',true,'GapOpen',100);
        A = A(1:2,:);
    %elseif (~isempty(seed))
    %    A = seedalign(kmers,seed);
    elseif (size(kmers,1)>2)
        A = multialign(kmers,'terminalGapAdjust',true,'GapOpen',100);
    else
        A = [];
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

% consensus
C = '';
for j = 1:size(N,2)
    [~,k] = max(N(:,j));
    C = [C L(k)];
end
C = [C{:}];

function A = seedalign(kmers,seed)

l = cellfun(@length,kmers);
f = cellfun(@min,strfind(kmers,seed));
s = max(f)-f;
e = l+s;
e = max(e)-e;

n = size(kmers,1);
for i = 1:n
    A(i,:) = [repmat('-',1,s(i)) kmers{i,:} repmat('-',1,e(i))];
end
