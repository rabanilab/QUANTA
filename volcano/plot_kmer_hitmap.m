function [h,I] = plot_kmer_hitmap(all_seq,E,P,fname,rowSort,thrP,thrE,mxP,mxE,I)

if (nargin < 5)
    rowSort = 0; % sort rows by: 0=effect-size; 1=pvalues; 2=kmers
end
if (nargin < 6)
    thrP = [];
end
if (nargin < 7)
    thrE = [];
end
if (nargin < 8)
    mxP = 30;
end
if (nargin < 9)
    mxE = 10;
end
if (nargin < 10)
    I = [];
end

% 0: red,black,green
% 1: red,white,blue
% 2: orange,white,purple
% 3: orange,white,blue
% 4: red,white,orange
% 5: brown,white,green
heatmap_color = 5;

% maximal number of rows to add text for kmers
maxText = 200;


% plot
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

[m,n] = size(E);

if (m>=3)
    N = align_kmers(all_seq);
    if (isempty(I))
        if (rowSort == 2)
            [~,i] = sortrows([E N*mxE/4]);
            all_seq = all_seq(i);
            P = P(i,:);
            E = E(i,:);
            N = N(i,:);
            [~,I] = cluster_sort([E N*mxE/4]);
        elseif (rowSort == 1)
            [~,i] = sortrows([E P]);
            all_seq = all_seq(i);
            P = P(i,:);
            E = E(i,:);
            N = N(i,:);
            [~,I] = cluster_sort([E P]);
        else
            [~,i] = sortrows(E);
            all_seq = all_seq(i);
            P = P(i,:);
            E = E(i,:);
            N = N(i,:);
            [~,I] = cluster_sort(E);
        end
    else
        i = (1:size(P,1))';
    end
    all_seq = all_seq(I);
    P = P(I,:);
    E = E(I,:);
    N = N(I,:);
    I = i(I);
elseif (m>0)
    N = align_kmers(all_seq);
else
    N = [];
end

ax = subplot(1,7,2:3);
imagesc(N);
colormap(ax,[1 1 1; 0 1 0; 0 0 1; 1 1 0; 1 0 0]);
c = colorbar('northOutside');
set(c,'ytick',0:4,'yticklabel',{'-' 'A' 'C' 'G' 'U'});
if (size(all_seq,1)<=maxText)
    set(gca,'ytick',1:size(all_seq,1),'yticklabel',all_seq);
else
    set(gca,'ytick',[]);
end
set(gca,'fontsize',14,'xtick',[]);
title('kmer sequence');

ax = subplot(1,7,4:5);
imagesc(E,[-1*mxE mxE]);
colormap(ax,gene_colormap(heatmap_color));
colorbar('northOutside');
if (isempty(thrE))
    title(sprintf('effect size (n=%d)',m));
else
    title(sprintf('effect size >= %.1f (n=%d)',thrE,m));
end
set(gca,'xtick',1:n,'xticklabel',regexprep(fname,'_','-'));
xtickangle(90);
set(gca,'ytick',[]);
set(gca,'fontsize',14);

ax = subplot(1,7,6:7);
imagesc(P,[-1*mxP mxP]);
colormap(ax,gene_colormap(heatmap_color));
colorbar('northOutside');
if (isempty(thrP))
    title(sprintf('-log10(p-value) (n=%d)',m));
else
    title(sprintf('-log10(p-value) >= %.1f (n=%d)',thrP,m));
end
set(gca,'xtick',1:n,'xticklabel',regexprep(fname,'_','-'));
xtickangle(90);
set(gca,'ytick',[]);
set(gca,'fontsize',14);



function [N,S] = align_kmers(all_seq)

n = size(all_seq,1);
if (n == 1)
    S = all_seq{1};
elseif (n == 2)
    S = multialign([all_seq;all_seq],'terminalGapAdjust',true,'GapOpen',100);
    S = S(1:2,:);
else
    S = multialign(all_seq,'terminalGapAdjust',true,'GapOpen',100);
end

N = zeros(size(S));
N(S=='A') = 1;
N(S=='C') = 2;
N(S=='G') = 3;
N(S=='T') = 4;
