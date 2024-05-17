function [W1,W2] = volcano_test(gid,rates,kname,Krange,ESIZE,ALPHA,Kdir,e_lim,p_lim,test,MHcorrect)
% sequence kmer analysis:
%   pvalue = difference of rate distribution between instances that contain each kmer and those that do not
%   effect = fold difference in median rates (log) between instances that contain each kmer and those that do not

if (nargin < 4)
    Krange = 3:7;
end
if (nargin < 5)
    ESIZE = 0.3;
end
if (nargin < 6)
    ALPHA = 0.01;
end
if (nargin < 8)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 7)
    e_lim = [-1 1];
end
if (nargin < 9)
    p_lim = [0 30];
end
if (nargin < 10)
    test = 'ks';
end
if (nargin < 11)
    MHcorrect = 'fdr';
end

% calculate pvalues
Yrows = [];
Q1 = [];
Q2 = [];
C = [];
E = [];
for k = Krange
    [e,p1,p2,c1,Xrows] = kmer_pvalue(k, gid, rates, kname, Kdir, test);
    Yrows = [Yrows; Xrows];
    Q1 = [Q1 p1];
    Q2 = [Q2 p2];
    C = [C c1];
    E = [E e];
end

% multiple hypothesis corrections
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

% significant p-values
i = (Q1<ALPHA);
write_text_file([kname '.' test '.smaller.txt'],sortrows(W(i,:),2));
W1 = W(i,:);
i = (Q1<ALPHA).*(E<=-1*ESIZE)==1;
write_text_file([kname '.' test '.e.smaller.txt'],sortrows(W(i,:),2));

i = (Q2<ALPHA);
write_text_file([kname '.' test '.larger.txt'],sortrows(W(i,:),3));
W2 = W(i,:);
i = (Q2<ALPHA).*(E>=ESIZE)==1;
write_text_file([kname '.' test '.e.larger.txt'],sortrows(W(i,:),3));

h = plot_volcano(E,Q1,Q2,Yrows,ALPHA,ESIZE,50,e_lim,p_lim);
xlabel('effect size (SMD)');
saveas(h, [kname '.' test '.all_seq.jpg'],'jpg');
saveas(h, [kname '.' test '.all_seq.eps'],'epsc');

close all;


function [e,p1,p2,c,Xrows] = kmer_pvalue(k, gid, rates, kname, SEQ_COUNT_DIR, test)

load([SEQ_COUNT_DIR '/kmers_counts.' num2str(k) '.mat'],'X','Xcols','Xrows');
Xcols = regexprep(Xcols,'.*;','');
fprintf('loaded data for k=%d: %d ids, %d kmers\n',k,size(X,2),size(X,1));

[~,j1,j2] = intersect(Xcols,gid);
Xcols = Xcols(j1);
X = X(:,j1);
gid = gid(j2(j2>0),:);
rates = rates(j2(j2>0),:);

i1 = full(sum(X,2)) > 0;
X = X(i1,:);
Xrows = Xrows(i1);

i2 = full(sum(X,1)) > 0;
X = X(:,i2);
gid = gid(i2);
rates = rates(i2);

fprintf('selected data for k=%d: %d ids, %d kmers\n',k,size(X,2),size(X,1));
min_test = 10;

if (exist([kname '.' num2str(k) test '.mat'],'file'))
    load([kname '.' num2str(k) test '.mat'],'e','p1','p2','c');
else
    % effect size
    e = [];
    sd = std(rates(~isnan(rates)),1);
    parfor i = 1:size(X,1)
        j = full((X(i,:)>0)');
        i0 = (j==0).*(~isnan(rates))==1;
        i1 = (j==1).*(~isnan(rates))==1;
        if ((sum(i0)>=min_test)*(sum(i1)>=min_test) == 1)
            e(i) = (mean(rates(i1)) - mean(rates(i0)))./sd;
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
        i0 = (j==0).*(~isnan(rates))==1;
        i1 = (j==1).*(~isnan(rates))==1;
        if ((sum(i0)>=min_test)*(sum(i1)>=min_test) == 1)
            c(i) = sum(i1);
            if (strcmp(test, 'mw')>0)
                p1(i) = ranksum(rates(i0)+E,rates(i1),'tail','right');
                p2(i) = ranksum(rates(i0)+E,rates(i1),'tail','left');
            elseif (strcmp(test, 'ks')>0)
                [~,p1(i)] = kstest2(rates(i0),rates(i1),'tail','smaller');
                [~,p2(i)] = kstest2(rates(i0),rates(i1),'tail','larger');
            end
        else
            c(i) = NaN;
            p1(i) = NaN;
            p2(i) = NaN;
        end
    end
    save([kname '.' num2str(k) test '.mat'],'e','p1','p2','c');
end
num2cell([k sum(isnan(e)) nanmean(round(e,2))])
