function [W1,W2,Ethr,W] = volcano_hyg(gid,labels,kname,Krange,ESIZE,ALPHA,Kdir,e_lim,p_lim,MHcorrect)
% sequence kmer analysis:
%   pvalue = difference of rate distribution between instances that contain each kmer and those that do not
%   effect = fold difference in median rates (log) between instances that contain each kmer and those that do not

if (nargin < 4)
    Krange = 3:7;
end
if (nargin < 5)
    ESIZE = 5;%0.3;
end
if (nargin < 6)
    ALPHA = 0.01;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    e_lim = [-1 1];
end
if (nargin < 9)
    p_lim = [0 30];
end
if (nargin < 10)
    MHcorrect = 'fdr';
end

% calculate pvalues
Yrows = [];
Q1 = [];
Q2 = [];
C = [];
E = [];
for k = Krange
    [e,p1,p2,c1,Xrows] = kmer_pvalue(k, gid, labels, kname, Kdir);
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

% effect size top percent
Ethr = prctile(abs(E),100-ESIZE);

% significant p-values
W = [Yrows num2cell([Q1;Q2;C;E])'];

i = (Q1<ALPHA);
write_text_file([kname '.hyg.smaller.txt'],sortrows(W(i,:),2));
W1 = W(i,:);
i = (Q1<ALPHA).*(E<=-1*Ethr)==1;
write_text_file([kname '.hyg.e.smaller.txt'],sortrows(W(i,:),2));

i = (Q2<ALPHA);
write_text_file([kname '.hyg.larger.txt'],sortrows(W(i,:),3));
W2 = W(i,:);
i = (Q2<ALPHA).*(E>=Ethr)==1;
write_text_file([kname '.hyg.e.larger.txt'],sortrows(W(i,:),3));

h = plot_volcano(E,Q1,Q2,Yrows,ALPHA,Ethr,0,e_lim,p_lim);
xlabel('effect size (%)');
axis square;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 10]); % x_width=10cm y_width=10cm
print([kname '.hyg.all.png'],'-dpng','-r300');
%saveas(h, [kname '.hyg.all.jpg'],'jpg');
h = plot_volcano(E,Q1,Q2,Yrows,ALPHA,Ethr,50,e_lim,p_lim);
xlabel('effect size (%)');
saveas(h, [kname '.hyg.all_seq.jpg'],'jpg');
%saveas(h, [kname '.hyg.all_seq.svg'],'svg');
close all;


function [e,p1,p2,c,Xrows] = kmer_pvalue(k, gid, labels, kname, SEQ_COUNT_DIR)

load([SEQ_COUNT_DIR '/kmers_counts.' num2str(k) '.mat'],'X','Xcols','Xrows');
Xcols = regexprep(Xcols,'.*;','');
fprintf('loaded data for k=%d: %d ids, %d kmers\n',k,size(X,2),size(X,1));

[~,j1,j2] = intersect(Xcols,gid);
Xcols = Xcols(j1);
X = X(:,j1);
gid = gid(j2(j2>0),:);
labels = labels(j2(j2>0),:);

i1 = full(sum(X,2)) > 0;
X = X(i1,:);
Xrows = Xrows(i1);

i2 = full(sum(X,1)) > 0;
X = X(:,i2);
gid = gid(i2);
labels = labels(i2);

fprintf('selected data for k=%d: %d ids, %d kmers\n',k,size(X,2),size(X,1));
min_test = 10;

if (exist([kname '.' num2str(k) 'hyg.mat'],'file'))
    load([kname '.' num2str(k) 'hyg.mat'],'e','p1','p2','c');
else
    % effect size
    e = [];
    parfor i = 1:size(X,1)
        j = full((X(i,:)>0)');
        i0 = (~isnan(labels));
        i1 = (j==1).*(~isnan(labels))==1;
        if ((sum(i0)>=min_test)*(sum(i1)>=min_test) == 1)
            e(i) = 100*sum(labels(i1))./sum(i1) - 100*sum(labels(i0))./sum(i0);
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
        i0 = (~isnan(labels));
        i1 = (j==1).*(~isnan(labels))==1;
        if ((sum(i0)>=min_test)*(sum(i1)>=min_test) == 1)
            c(i) = sum(i1);
            p1(i) = hypergeometric_pvalue(sum(labels(i1)),sum(i0),sum(i1),sum(labels(i0)),1);
            p2(i) = hypergeometric_pvalue(sum(labels(i1)),sum(i0),sum(i1),sum(labels(i0)),0);
        else
            c(i) = NaN;
            p1(i) = NaN;
            p2(i) = NaN;
        end
    end

    % save results
    save([kname '.' num2str(k) 'hyg.mat'],'e','p1','p2','c');
end
num2cell([k sum(isnan(e)) nanmean(round(e,2))])
