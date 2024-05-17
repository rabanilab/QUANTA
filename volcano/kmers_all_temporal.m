function kmers_all_temporal(Xids,Dname,Data,output_dir,KMER_RANGE,maxK,Esize,alpha,Erange,Prange,Kdir,findK)

if (nargin < 5)
    KMER_RANGE = 3:7;
end
if (nargin < 6)
    maxK = 200;
end
if (nargin < 7)
    Esize = 2;
end
if (nargin < 8)
    alpha = 0.01;
end
if (nargin < 9)
    Erange = [-10 10];
end
if (nargin < 10)
    Prange = [0 30];
end
if (nargin < 11)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 12)
    findK = 0;
end

n = size(Data,2);
mkdir(output_dir);

% KS test to identify significant kmers
if (exist([output_dir '/ks.mat'],'file'))
    load([output_dir '/ks.mat'],'W1','W2');
    fprintf('load KS test data\n');
else
    W1 = cell(1,n);
    W2 = cell(1,n);
    for i = 1:n
        Dname(i)
        I = isnan(Data(:,i)) == 0;
        [W1{i},W2{i}] = volcano_ks(Xids(I),Data(I,i),[output_dir '/' Dname{i}],KMER_RANGE,Esize,alpha,Erange,Kdir,Prange);
    end
    save([output_dir '/ks.mat'],'W1','W2');
end

% select thresholds
[~,~,~,thrP1,thrE1,mxP1,mxE1] = kmer_merge_temporal(W1,2,alpha,Esize,maxK);
[~,~,~,thrP2,thrE2,mxP2,mxE2] = kmer_merge_temporal(W2,3,alpha,Esize,maxK);
thrE = max(thrE1,thrE2);
thrP = max(thrP1,thrP2);
mxP = max(mxP1,mxP2);
if (mxP <= 0)
    mxP = 0.1;
end
mxE = max(mxE1,mxE2);
if (mxE <= 0)
    mxE = 0.1;
end
num2cell([thrE thrP mxE mxP])

% plot matrices
[all_seq1,P1,E1] = kmer_merge_temporal(W1,2,alpha,Esize,maxK,thrP1,thrE1,mxP,mxE);
h = plot_kmer_hitmap(all_seq1,E1,P1,Dname,1,thrP,thrE,mxP,mxE);
saveas(h, [output_dir '/matrix.smaller.jpg'],'jpg');
write_text_file([output_dir '/matrix_P.smaller.txt'],[all_seq1 num2cell(P1)]);
write_text_file([output_dir '/matrix_E.smaller.txt'],[all_seq1 num2cell(E1)]);
close all;

[all_seq2,P2,E2] = kmer_merge_temporal(W2,3,alpha,Esize,maxK,thrP2,thrE2,mxP,mxE);
h = plot_kmer_hitmap(all_seq2,E2,P2,Dname,1,thrP,thrE,mxP,mxE);
saveas(h, [output_dir '/matrix.larger.jpg'],'jpg');
write_text_file([output_dir '/matrix_P.larger.txt'],[all_seq2 num2cell(P2)]);
write_text_file([output_dir '/matrix_E.larger.txt'],[all_seq2 num2cell(E2)]);
close all;

% find groups of k-mers
if (findK)
    pscore1 = max(E1,[],2);
    if (size(pscore1,1)>=5)
        [PWM1,Q1,PWMseq1] = merge_kmers_by_kmeans(all_seq1,pscore1);
        if (~isempty(PWM1))
            write_text_file([output_dir '/pwm.smaller.txt'], sortrows(PWMseq1,[4 5]));
        end
        for i = 1:size(PWM1,1)
            h = pwm_logo(PWM1{i});
            saveas(h,[output_dir '/pwm.smaller.' Q1{i,1} '.' num2str(Q1{i,2}) '.jpg'],'jpg');
            close(h);
        end
    end
    close all;
    
    pscore2 = max(E2,[],2);
    if (size(pscore2,1)>=5)
        [PWM2,Q2,PWMseq2] = merge_kmers_by_kmeans(all_seq2,pscore2);
        if (~isempty(PWM2))
            write_text_file([output_dir '/pwm.larger.txt'], sortrows(PWMseq2,[4 5]));
        end
        for i = 1:size(PWM2,1)
            h = pwm_logo(PWM2{i});
            saveas(h,[output_dir '/pwm.larger.' Q2{i,1} '.' num2str(Q2{i,2}) '.jpg'],'jpg');
            close(h);
        end
    end
    close all;
end

function [all_seq,P,E,thrP,thrE,mxP,mxE] = kmer_merge_temporal(W,k,alpha,Esize,maxK,thrP,thrE,mxP,mxE,dMAX)
% select by maximal pvalue

if (nargin < 5)
    maxK = 400;
end
if (nargin < 6)
    thrP = -1*log10(alpha);
end
if (nargin < 7)
    thrE = Esize;
end
if (nargin < 8)
    mxP = 30;
end
if (nargin < 9)
    mxE = 10;
end
if (nargin < 10)
    dMAX = 0.01;
end

% merge all kmers
n = size(W,2);
all_seq = [];
for i = 1:n
    if (~isempty(W{i}))
        all_seq = [all_seq;W{i}(:,1)];
    end
end
all_seq = unique(all_seq);

all_P = ones(size(all_seq,1),n);
all_E = zeros(size(all_seq,1),n);
for i = 1:n
    if (~isempty(W{i}))
        [w1,w2] = ismember(all_seq,W{i}(:,1)); % [all_seq(w1) W1{i}(w2(w2>0),1)]
        if (sum(w1)>0)
            all_P(w1,i) = cell2mat(W{i}(w2(w2>0),k));
            all_E(w1,i) = cell2mat(W{i}(w2(w2>0),5));
        end
    end
end
all_P(all_P<1e-100) = 1e-100;
E = abs(all_E);
P = -1*log10(all_P);

% select significant kmers
[allP,i] = max(P,[],2);
allE = get_max(E,i);

i = 1;
I = (allP>=thrP).*(allE>=thrE)==1;
while (sum(I)>maxK)
    thrP = min(allP(I));
    thrE = min(allE(I));
    %if (i/2 == round(i/2))
    %    thrE = thrE + dMAX;
    %else
        thrP = thrP + dMAX;
    %end
    i = i+1;
    I = (allP>=thrP).*(allE>=thrE)==1;
    fprintf('n=%d [%.2f, %.2f]\n',sum(I),thrP,thrE);
end
all_seq = all_seq(I);
P = P(I,:);
E = E(I,:);

% define thresholds
if (isempty(P))
    mx = 0;
else
    mx = prctile(P(:),95);
end
if (mxP > mx)
    mxP = mx;
end

if (isempty(E))
    mx = 0;
else
    mx = prctile(E(:),95);
end
if (mxE > mx)
    mxE = mx;
end
fprintf('select temporal: %d kmers in matrix [%.1f, %.1f]\n',sum(I),thrP,thrE);


function M = get_max(X,i)

X = X';
i = i';

n = size(X,1);
M = X(i+cumsum(n*ones(size(i)))-n);
M = M';
