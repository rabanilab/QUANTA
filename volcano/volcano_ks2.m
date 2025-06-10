function [W1, W2] = volcano_ks2(gid,rate1,rate2,kname,Krange,ESIZE,ALPHA,Kdir,e_lim,p_lim,MHcorrect)
% sequence kmer analysis:
%   pvalue = difference of rate distribution between instances that contain each kmer in two groups
%   effect = fold difference in median rates (log) between instances that contain each kmer in two groups

if (nargin < 5)
    Krange = 3:7;
end
if (nargin < 6)
    ESIZE = 0.3;
end
if (nargin < 7)
    ALPHA = 0.01;
end
if (nargin < 8)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 9)
    e_lim = [-1 1];
end
if (nargin < 10)
    p_lim = [0 30];
end
if (nargin < 11)
    MHcorrect = 'fdr';
end

Yrows = [];
Q1 = [];
Q2 = [];
C = [];
E = [];
for k = Krange
    [e,c,p1,p2,Xrows] = kmer_pvalue(k, gid, rate1, rate2, kname, Kdir);
    Yrows = [Yrows; Xrows];
    Q1 = [Q1 p1];
    Q2 = [Q2 p2];
    C = [C c];
    E = [E e];
end

% all kmers in range
if (strcmp(MHcorrect,'bf'))
    b = max(size(Q1(:))) + max(size(Q2(:)));
    Q1 = Q1*b;
    Q1(Q1>1) = 1;
    Q2 = Q2*b;
    Q2(Q2>1) = 1;
else
    Q1 = mafdr(Q1,'BHFDR',true);
    Q2 = mafdr(Q2,'BHFDR',true);
end
W = [Yrows num2cell([Q1;Q2;C;E])'];

i = (Q1<ALPHA);
write_text_file([kname '.ks2.larger.txt'],sortrows(W(i,:),2));
W1 = W(i,:);
i = (Q1<ALPHA).*(E>=ESIZE)==1;
write_text_file([kname '.ks2.e.larger.txt'],sortrows(W(i,:),2));

i = (Q2<ALPHA);
write_text_file([kname '.ks2.smaller.txt'],sortrows(W(i,:),3));
W2 = W(i,:);
i = (Q2<ALPHA).*(E<=-1*ESIZE)==1;
write_text_file([kname '.ks2.e.smaller.txt'],sortrows(W(i,:),3));

h = plot_volcano(E,Q1,Q2,Yrows,ALPHA,ESIZE,50,e_lim,p_lim);
xlabel('effect size (fold-difference #1 - #2, log2)');
title('all kmers');
saveas(h, [kname '.ks2.all_seq.jpg'],'jpg');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 10]); % x_width=10cm y_width=10cm
saveas(h, [kname '.ks2.all_seq.svg'],'svg');
%saveas(h, [kname '.ks2.all_seq.eps'],'epsc');

close all;


function [e,c,p1,p2,Xrows] = kmer_pvalue(k, gid, rate1, rate2, kname, SEQ_COUNT_DIR)

min_test = 10;

load([SEQ_COUNT_DIR '/kmers_counts.' num2str(k) '.mat'],'X','Xcols','Xrows');
Xcols = regexprep(Xcols,'.*;','');

[~,j1,j2] = intersect(Xcols,gid); % sum(strcmp(Xcols(j1),gid(j2(j2>0))))
Xcols = Xcols(j1);
X = X(:,j1);
i1 = full(sum(X,2))>0;
X = X(i1,:);
Xrows = Xrows(i1);
gid = gid(j2(j2>0),:);
rate1 = rate1(j2(j2>0),:);
rate2 = rate2(j2(j2>0),:);
fprintf('intersect k=%d: %d ids overlap, %d kmers\n',k,max(size(gid)),size(X,1));

if (exist([kname '.' num2str(k) 'ks2.mat'],'file'))
    load([kname '.' num2str(k) 'ks2.mat'],'e','p1','p2','c');
else
    % effect size
    e = [];
    sd = std([rate1(~isnan(rate1));rate2(~isnan(rate2))],1);
    parfor i = 1:size(X,1)
        j = full((X(i,:)>0)');
        i1 = (~isnan(rate1)).*(~isnan(rate2)).*(j==1) == 1;
        if (sum(i1)>=min_test)
            e(i) = (mean(rate1(i1)) - mean(rate2(i1)))./sd;
        else
            e(i) = NaN;
        end
    end
    
    % pvalue
    p1 = [];
    p2 = [];
    c = [];
    parfor i = 1:size(X,1)
        j = full((X(i,:)>0)');
        i1 = (~isnan(rate1)).*(~isnan(rate2)).*(j==1) == 1;
        if (sum(i1)>=min_test)
            c(i) = sum(i1);
            [~,p1(i)] = kstest2(rate1(i1),rate2(i1),'tail','smaller');
            [~,p2(i)] = kstest2(rate1(i1),rate2(i1),'tail','larger');
        else
            c(i) = NaN;
            p1(i) = NaN;
            p2(i) = NaN;
        end
    end
    save([kname '.' num2str(k) 'ks2.mat'],'e','p1','p2','c');
end

num2cell([k sum(isnan(e)) nanmean(round(e,2))])
