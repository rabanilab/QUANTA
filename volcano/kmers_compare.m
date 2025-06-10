function kmers_compare(Dids,Dseq,Dval,Dnames,output_dir,Krange,Kdir,esize,palpha,erange,prange,norm_len)

if (nargin < 6)
    Krange = 3:7;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    esize = 5;%0.3;
end
if (nargin < 9)
    palpha = 0.01;
end
if (nargin < 10)
    erange = [-1 1];
end
if (nargin < 11)
    prange = [0 10];
end
if (nargin < 13)
    norm_len = 1;
end

n = size(Dval,2);
Dlen = cellfun(@length,Dseq);
len_fold = 10;

mkdir(output_dir);


% plot data
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
mnD = prctile(Dval(:),1);
mxD = prctile(Dval(:),99);
x = mnD:((mxD-mnD)/25):mxD;
y = hist(Dval,x);
plot(x,y./sum(y),'-','linewidth',2);
axis tight;
xlabel('input values');
ylabel('fraction');
set(gca,'fontsize',18);
L = regexprep(strcat(Dnames',' (n=',num2str(sum(~isnan(Dval))'),')'),'n=  *','n=');
legend(L,'location','bestOutside','box','off');
saveas(h, [output_dir '/data.jpg'],'jpg');

close all;

% ---------------------------------------------------------------------
% KS test to identify significant kmers
% W1/W2 = [kmer] [pvalue neg. E] [pvalue pos. E] [gene +kmer count] [effect size]
% ---------------------------------------------------------------------
W = cell(1,n);
Ethr = zeros(1,n);
for i = 1:n
    % remove extreme lengths: very long/short sequences skew the statistical test
    I = (~isnan(Dval(:,i))).*(Dlen > (1/len_fold)*median(Dlen)).*(Dlen < len_fold*median(Dlen)) == 1;

    % test if there is a correlation between length and property
    p = fit_linear_model(Dlen(I),Dval(I,i));
    x = [1 max(Dlen(I))];
    y = polyval(p,x);
    zy = y - (polyval(p,x)-p(2));
    zDval = Dval(I,i) - (polyval(p,Dlen(I)) - p(2));

    h = figure;
    subplot(1,2,1);
    hold on;
    w = 0.01*randn(size(Dlen(I)));
    dscatter(Dlen(I)+w,Dval(I,i),'MSIZE',25);
    plot(x,y,'-k');
    hold off;
    axis tight;
    axis square;
    xlabel('sequence length (nt)');
    ylabel(Dnames{i});
    c1 = corr(Dlen(I),Dval(I,i));
    title(sprintf('n=%d, c=%.2f',sum(I),c1));
    subplot(1,2,2);
    hold on;
    dscatter(Dlen(I)+w,zDval,'MSIZE',25);
    plot(x,zy,'-k');
    hold off;
    axis tight;
    axis square;
    xlabel('sequence length (nt)');
    ylabel(Dnames{i});
    c2 = corr(Dlen(I),zDval);
    title(sprintf('n=%d, c=%.2f',sum(I),c2));
    saveas(h,[output_dir '/corr.' Dnames{i} '.jpg'],'jpg');
    close all;

    % KS test
    if (norm_len)
        [~,~,Ethr(i),W{i}] = volcano_ks(Dids(I),zDval,[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
    else
        [~,~,Ethr(i),W{i}] = volcano_ks(Dids(I),Dval(I,i),[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
    end
end

% merge kmer enrichments
Kall = [];
for i = 1:n
    Kall = [Kall;W{i}(:,1)];
end
Kall = unique(Kall);
m = size(Kall,1);

Wp = ones(m,n);
We = zeros(m,n);
for i = 1:n
    e = cell2mat(W{i}(:,5)) > 0;
    [~,j] = ismember(W{i}(e,1),Kall);
    Wp(j,i) = cell2mat(W{i}(e,3));
    We(j,i) = cell2mat(W{i}(e,5));
    e = cell2mat(W{i}(:,5)) < 0;
    [~,j] = ismember(W{i}(e,1),Kall);
    Wp(j,i) = -1*cell2mat(W{i}(e,2));
    We(j,i) = cell2mat(W{i}(e,5));
end
j1 = sum((abs(Wp)<palpha).*(abs(We)<max(Ethr)) == 1,2) > 0;
Wp = Wp(j1,:);
logWp = -1*sign(Wp).*log10(abs(Wp));
We = We(j1,:);
K1 = Kall(j1);

write_text_file([output_dir '/ks.merge.P.txt'],[['kmer' Dnames];[K1 num2cell(logWp)]]);
write_text_file([output_dir '/ks.merge.E.txt'],[['kmer' Dnames];[K1 num2cell(We)]]);


% plot all comparisons
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
nr = ceil(n/3);

% pvalues
for i1 = 1:n
    clf;
    for i2 = 1:n
        subplot(3,nr,i2);
        j = max(abs(logWp(:,[i1 i2])),[],2) > -1*log10(palpha);
        mn = min([logWp(j,i1);logWp(j,i2)]);
        mx = max([logWp(j,i1);logWp(j,i2)]);
        hold on;
        dscatter(logWp(j,i1),logWp(j,i2),'MSIZE',25);
        line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
        line([0 0],[mn mx],'LineStyle','-','color','k','linewidth',1);
        line([mn mx],[0 0],'LineStyle','-','color','k','linewidth',1);
        hold off;
        axis square;
        axis tight;
        box on;
        xlabel(Dnames{i1});
        ylabel(Dnames{i2});
        title(sprintf('n=%d,r=%.2f [%.1f,%.1f]',sum(j),corr(logWp(j,i1),logWp(j,i2)),mn,mx));
        set(gca,'fontsize',14,'xtick',[],'ytick',[]);
    end
    saveas(h, [output_dir '/kmer_corr.' Dnames{i1} '.p.jpg'],'jpg');
end

% effect size
for i1 = 1:n
    clf;
    for i2 = 1:n
        subplot(3,nr,i2);
        j = max(abs(logWp(:,[i1 i2])),[],2) > -1*log10(palpha);
        mn = min([We(j,i1);We(j,i2)]);
        mx = max([We(j,i1);We(j,i2)]);
        hold on;
        dscatter(We(j,i1),We(j,i2),'MSIZE',25);
        line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
        line([mn mx],1.5*[mn mx],'LineStyle','-','color','r','linewidth',1);
        line(1.5*[mn mx],[mn mx],'LineStyle','-','color','r','linewidth',1);
        line([0 0],[mn mx],'LineStyle','-','color','k','linewidth',1);
        line([mn mx],[0 0],'LineStyle','-','color','k','linewidth',1);
        hold off;
        axis square;
        axis tight;
        box on;
        xlabel(Dnames{i1});
        ylabel(Dnames{i2});
        title(sprintf('n=%d,r=%.2f [%.1f,%.1f]',sum(j),corr(We(j,i1),We(j,i2)),mn,mx));
        set(gca,'fontsize',14,'xtick',[],'ytick',[],'ylim',[mn mx],'xlim',[mn mx]);
    end
    saveas(h, [output_dir '/kmer_corr.' Dnames{i1} '.e.jpg'],'jpg');
end

% plot pairwise comparisons
n2 = ceil(n/2);
nr = ceil(n2/3);

clf;
for i1 = 1:n2
    subplot(3,nr,i1);
    i2 = i1+n2;
    j = max(abs(logWp(:,[i1 i2])),[],2) > -1*log10(palpha);
    mn = min([logWp(j,i1);logWp(j,i2)]);
    mx = max([logWp(j,i1);logWp(j,i2)]);
    hold on;
    dscatter(logWp(j,i1),logWp(j,i2),'MSIZE',25);
    line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
    line([0 0],[mn mx],'LineStyle','-','color','k','linewidth',1);
    line([mn mx],[0 0],'LineStyle','-','color','k','linewidth',1);
    hold off;
    axis square;
    axis tight;
    box on;
    xlabel(Dnames{i1});
    ylabel(Dnames{i2});
    title(sprintf('n=%d,r=%.2f [%.1f,%.1f]',sum(j),corr(logWp(j,i1),logWp(j,i2)),mn,mx));
    set(gca,'fontsize',14,'xtick',[],'ytick',[]);
end
saveas(h, [output_dir '/kmer_corr.pw.jpg'],'jpg');

clf;
L = {'kmer' 'type1' 'esize1' 'type2' 'esize2'};
for i1 = 1:n2
    subplot(3,nr,i1);
    i2 = i1+n2;
    j = max(abs(logWp(:,[i1 i2])),[],2) > -1*log10(palpha);
    mn = min([We(j,i1);We(j,i2)]);
    mx = max([We(j,i1);We(j,i2)]);
    hold on;
    dscatter(We(j,i1),We(j,i2),'MSIZE',25);
    line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
    line([mn mx],1.5*[mn mx],'LineStyle','-','color','r','linewidth',1);
    line(1.5*[mn mx],[mn mx],'LineStyle','-','color','r','linewidth',1);
    line([0 0],[mn mx],'LineStyle','-','color','k','linewidth',1);
    line([mn mx],[0 0],'LineStyle','-','color','k','linewidth',1);
    hold off;
    axis square;
    axis tight;
    box on;
    xlabel(Dnames{i1});
    ylabel(Dnames{i2});
    title(sprintf('n=%d,r=%.2f [%.1f,%.1f]',sum(j),corr(We(j,i1),We(j,i2)),mn,mx));
    set(gca,'fontsize',14,'xtick',[],'ytick',[],'ylim',[mn mx],'xlim',[mn mx]);
    k = j.*((abs(We(:,i1))>=1.5*abs(We(:,i2)))+(abs(We(:,i2))>=1.5*abs(We(:,i1)))+(sign(We(:,i2))~=sign(We(:,i1))) > 0) == 1;
    %plot(We(k,i1),We(k,i2),'.r');
    L = [L;[K1(k) repmat(Dnames(i1),sum(k),1) num2cell(We(k,i1)) repmat(Dnames(i2),sum(k),1) num2cell(We(k,i2))]];
end
saveas(h, [output_dir '/kmer_corr.ew.jpg'],'jpg');
write_text_file([output_dir '/kmer_corr.ew.txt'],L);

close all;


function p = fit_linear_model(x,y)

u1 = unique(x);
if (max(size(u1)) <= 2)
    p = [0 mean(y)];
else
    %p = polyfit(x,y,1);
    p = robustfit(x,y);
    p = p(end:-1:1)';
end

