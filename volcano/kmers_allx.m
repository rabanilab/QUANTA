function kmers_allx(Dids,Dseq,Dval,Dnames,output_dir,Krange,Kdir,esize,palpha,erange,prange,sort_type)

if (nargin < 6)
    Krange = 3:7;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    esize = 0.5;
end
if (nargin < 9)
    palpha = 0.01;
end
if (nargin < 10)
    erange = [-1 1];
end
if (nargin < 11)
    prange = [0 30];
end
if (nargin < 12)
    sort_type = 1; % 1 = by pvalues; 2 = by kmers
end

maxD = 3;
minP = 1e-200;
min_overlap = min(Krange);
n = size(Dval,2);

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
legend(Dnames,'location','bestOutside','box','off');
saveas(h, [output_dir '/data.jpg'],'jpg');

fprintf('nan values:\n');
[Dnames' num2cell(sum(isnan(Dval))')]

% KS test to identify significant kmers
if (exist([output_dir '/ks.mat'],'file'))
    load([output_dir '/ks.mat'],'W1','W2');
    fprintf('load KS test data\n');
else
    W1 = cell(1,n);
    W2 = cell(1,n);
    for i = 1:n
        Dnames{i}
        I = isnan(Dval(:,i)) == 0;
        if (sum(I)>10)
            [W1{i},W2{i}] = volcano_ks(Dids(I),Dval(I,i),[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
        else
            W1{i} = [];
            W2{i} = [];
        end
    end
    save([output_dir '/ks.mat'],'W1','W2');
end

% for each experiment, merge kmers by sequence overlaps
% [1] collect kmers with sequence overlaps, combine into a PWM
% [2] score all input sequences by PWM, exclude PWMs with non-significant correlation to values
if (exist([output_dir '/pwm.mat'],'file'))
    load([output_dir '/pwm.mat'],'P1','P2');
    fprintf('load PWM data\n');
else
    S1 = cell(1,n);
    P1 = cell(1,n);
    C1 = cell(1,n);
    S2 = cell(1,n);
    P2 = cell(1,n);
    C2 = cell(1,n);
    for j = 1:n
        Dnames{j}
        if (~isempty(W1{j}))
            [S1{j},P1{j}] = merge_kmers_by_seq(W1{j}(:,1),cell2mat(W1{j}(:,2)),maxD,min_overlap);
            [P1{j},S1{j},C1{j}] = pwm_filter_corr(P1{j},S1{j},Dseq,Dval(:,j),palpha,'left');
        end
        if (~isempty(W2{j}))
            [S2{j},P2{j}] = merge_kmers_by_seq(W2{j}(:,1),cell2mat(W2{j}(:,3)),maxD,min_overlap);
            [P2{j},S2{j},C2{j}] = pwm_filter_corr(P2{j},S2{j},Dseq,Dval(:,j),palpha,'right');
        end
        h = pwm_plot_corr(C1{j},P1{j},S1{j},C2{j},P2{j},S2{j});
        saveas(h, [output_dir '/pwm.' Dnames{j} '.1.jpg'],'jpg');
        close all;
    end
    save([output_dir '/pwm.mat'],'P1','P2','S1','S2','C1','C2');
end

% combine PWMs from all different datasets:
% [1] combine all datasets 
% [2] exclude similar PWMs
% [3] score all sequences by PWMs
if (exist([output_dir '/pwm_all.mat'],'file'))
    load([output_dir '/pwm_all.mat'],'S_all','E_all','Q_all');
    fprintf('load PWM combined\n');
else
    P_all = [];
    for i = 1:max(size(P1))
        P_all = [P_all P1{i}];
    end
    for i = 1:max(size(P2))
        P_all = [P_all P2{i}];
    end
    
    if (~isempty(P_all))
        [P_all,S_all] = pwm_filter_similar(P_all,maxD,min_overlap);
        C_all = pwm_score_all(P_all,Dseq);
        [E_all,Q_all] = corr(Dval,C_all,'Rows','pairwise');
        q = mafdr(Q_all(:),'BHFDR',true);
        Q_all = reshape(q,size(Q_all))';
        E_all = E_all';
        E_all(isnan(E_all)) = 0;
    else
        S_all = [];
        C_all = [];
        Q_all = [];
        E_all = [];
    end
    save([output_dir '/pwm_all.mat'],'S_all','E_all','Q_all','C_all','P_all');
end

if (isempty(Q_all))
    return;
end

% plot results
P_all = -1*log10(Q_all);
P_all(P_all > -1*log10(minP)) = -1*log10(minP);
P_all(isnan(P_all)) = 0;
mxP = max(P_all(:));
mxE = max(abs(E_all(:)));

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) 0.5*scrsz(3) scrsz(4)]);
for j = 1:n
    clf;
    hold on;
    for i = 1:size(E_all,1)
        if (P_all(i,j) >= -1*log10(palpha))
            seq = regexprep(S_all(i),'T','U');
            text(E_all(i,j),P_all(i,j),seq,'fontsize',14);
        end
    end
    n1 = 0.01*randn(size(E_all(:,j)));
    n2 = 0.01*randn(size(P_all(:,j)));
    dscatter(E_all(:,j)+n1,P_all(:,j)+n2,'MSIZE',50);
    line([-1*mxE mxE],-1*log10([palpha palpha]),'LineStyle','-','color','b','linewidth',1);
    line([0 0],[0 mxP],'LineStyle','-','color','k','linewidth',1);
    hold off;
    xlabel('effect size');
    ylabel('p-value');
    title(Dnames{j});
    set(gca,'xlim',[-1*mxE mxE],'ylim',[0 mxP]);
    saveas(h, [output_dir '/matrix.' Dnames{j} '.volcano.jpg'],'jpg');
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
Cp = ones(n,n);
Ce = ones(n,n);
for j = 1:n
    clf;
    w = 1;
    o = P_all(:,j) >= -1*log10(palpha);
    if (sum(o) > 0)
        for i = setdiff(1:n,j)
            Cp(i,j) = corr(P_all(o,j),P_all(o,i));
            Ce(i,j) = corr(E_all(o,j),E_all(o,i));
            subplot(2,n-1,w);
            hold on;
            %for k = find(o)'
            %    text(P_all(k,j),P_all(k,i),S_all(k,1),'fontsize',10);
            %end
            if (sum(o) > 5)
                t = 0.001*randn(sum(o),1);
                dscatter(P_all(o,j)+t,P_all(o,i)+t,'MSIZE',50);
            else
                plot(P_all(o,j),P_all(o,i),'.','markersize',20);
            end
            line([0 mxP],[0 mxP],'LineStyle','-','color','k','linewidth',1);
            line([0 mxP],-1*log10([palpha palpha]),'LineStyle','-','color','b','linewidth',1);
            line(-1*log10([palpha palpha]),[0 mxP],'LineStyle','-','color','b','linewidth',1);
            hold off;
            box on;
            xlabel([Dnames{j} ' (p-value)']);
            ylabel([Dnames{i} ' (p-value)']);
            axis square;
            set(gca,'xlim',[0 mxP],'ylim',[0 mxP]);
            title(sprintf('n=%d, r=%.2f',sum(o),Cp(i,j)));
            subplot(2,n-1,n-1+w);
            hold on;
            %for k = find(o)'
            %    text(E_all(k,j),E_all(k,i),S_all(k,1),'fontsize',10);
            %end
            if (sum(o) > 5)
                t = 0.001*randn(sum(o),1);
                dscatter(E_all(o,j)+t,E_all(o,i)+t,'MSIZE',50);
            else
                plot(E_all(o,j),E_all(o,i),'.','markersize',20);
            end
            line([-1*mxE mxE],[-1*mxE mxE],'LineStyle','-','color','k','linewidth',1);
            line([0 0],[-1*mxE mxE],'LineStyle','-','color','k','linewidth',1);
            line([-1*mxE mxE],[0 0],'LineStyle','-','color','k','linewidth',1);
            hold off;
            box on;
            xlabel([Dnames{j} ' (effect)']);
            ylabel([Dnames{i} ' (effect)']);
            axis square;
            set(gca,'xlim',[-1*mxE mxE],'ylim',[-1*mxE mxE]);
            title(sprintf('n=%d, r=%.2f',sum(o),Ce(i,j)));
            w = w + 1;
        end
    end
    saveas(h, [output_dir '/matrix.' Dnames{j} '.corr.jpg'],'jpg');
end
Cp(isnan(Cp)) = 0;
Ce(isnan(Ce)) = 0;

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
[~,i] = cluster_sort(Ce);
subplot(1,2,1);
mx = max(abs(Cp(:)));
imagesc(Cp(i,i),[-1*mx mx]);
axis square;
colormap(gene_colormap(1));
colorbar;
set(gca,'xtick',1:max(size(Dnames)),'xticklabels',Dnames(i),'fontsize',16);
xtickangle(90);
set(gca,'ytick',1:max(size(Dnames)),'yticklabels',Dnames(i));
ylabel('comparison base');
title('correlation of p-values between motifs');
subplot(1,2,2);
mx = max(abs(Ce(:)));
imagesc(Ce(i,i),[-1*mx mx]);
axis square;
colormap(gene_colormap(1));
colorbar;
set(gca,'xtick',1:max(size(Dnames)),'xticklabels',Dnames(i),'fontsize',16);
xtickangle(90);
set(gca,'ytick',1:max(size(Dnames)),'yticklabels',Dnames(i));
ylabel('comparison base');
title('correlation of effect sizes between motifs');
saveas(h, [output_dir '/matrix.corr.jpg'],'jpg');

mxP = 10;%prctile(P_all(:),98);30;
mxE = 0.05;%prctile(abs(E_all(:)),98);0.1;
P = P_all;
P(E_all<0) = -1*P_all(E_all<0);
[h,i] = plot_kmer_hitmap(S_all,E_all,P,Dnames,sort_type,-1*log10(palpha),esize,mxP,mxE);
saveas(h, [output_dir '/matrix.all.jpg'],'jpg');
write_text_file([output_dir '/matrix.P.all.txt'],[S_all(i) num2cell(P(i,:))]);
write_text_file([output_dir '/matrix.E.all.txt'],[S_all(i) num2cell(E_all(i,:))]);

k = zeros(size(P,1),1);
for j = 1:n
    pj = abs(P(:,j));
    k = k + (pj > prctile(pj,90)) > 0;
end
Sk = S_all(k,:);
Ek = E_all(k,:);
Pk = P(k,:);
[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),esize,mxP,mxE);
saveas(h, [output_dir '/matrix.90.jpg'],'jpg');
write_text_file([output_dir '/matrix.P.90.txt'],[Sk(i) num2cell(Pk(i,:))]);
write_text_file([output_dir '/matrix.E.90.txt'],[Sk(i) num2cell(Ek(i,:))]);

for j = 1:n
    pj = abs(P(:,j));
    k = (pj > prctile(pj,90)) > 0;
    Sk = S_all(k,:);
    Ek = E_all(k,:);
    Pk = P(k,:);
    [~,I] = sort(Pk(:,j));
    [h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),esize,mxP,mxE,I);
    saveas(h, [output_dir '/matrix.' Dnames{j} '.jpg'],'jpg');
    write_text_file([output_dir '/matrix.P.' Dnames{j} '.txt'],[Sk(i) num2cell(Pk(i,:))]);
    write_text_file([output_dir '/matrix.E.' Dnames{j} '.txt'],[Sk(i) num2cell(Ek(i,:))]);
end

close all;


function [P1,S1,SCORE,POS] = pwm_filter_corr(P1,S1,Dseq,Dval,alpha,tail_type)

fprintf('Input: %d PWMs\n',size(S1,1));

[SCORE,POS] = pwm_score_all(P1,Dseq);
[~,Q1] = corr(Dval,SCORE,'Rows','pairwise','tail',tail_type);
q = mafdr(Q1(:),'BHFDR',true);
Q1 = reshape(q,size(Q1))';
Q1 = min(Q1,[],2);
S1 = [S1 num2cell(Q1)];

I = Q1<alpha;
fprintf('Final: %d PWMs\n',sum(I));
%S1(I==0,:)

S1 = S1(I,:);
P1 = P1(I);
SCORE = SCORE(:,I);
POS = POS(:,I);

function [h,i1,i2] = pwm_plot_corr(C1,P1,S1,C2,P2,S2)

if (~isempty(C1))
    CC1 = corr(C1);
else
    CC1 = [];
    i1 = [];
end
if (~isempty(C2))
    CC2 = corr(C2);
else
    CC2 = [];
    i2 = [];
end

p1 = prctile([CC1(:);CC2(:)],95);
p2 = prctile([CC1(:);CC2(:)],5);
mnC = min(p1,p2);
mxC = max(p1,p2);
if (mnC == mxC)
    mxC = mxC+0.1;
end
if ((mxC>0)&&(mnC<0))
    mnC = -1*max(abs(mnC),mxC);
    mxC = max(abs(mnC),mxC);
    c = gene_colormap(1);
else
    c = intensity_colormap();
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(1,2,1);
if (~isempty(C1))
    [~,i1] = cluster_sort(CC1);
    imagesc(CC1(i1,i1),[mnC mxC]);
    set(gca,'ytick',1:size(C1,2),'yticklabel',S1(i1,1));
    set(gca,'xtick',[]);
end
axis square;
colormap(c);
colorbar;
title(sprintf('smaller (n=%d)',max(size(P1))));
subplot(1,2,2);
if (~isempty(C2))
    [~,i2] = cluster_sort(CC2);
    imagesc(CC2(i2,i2),[mnC mxC]);
    set(gca,'ytick',1:size(C2,2),'yticklabel',S2(i2,1));
    set(gca,'xtick',[]);
end
axis square;
colormap(c);
colorbar;
title(sprintf('larger (n=%d)',max(size(P2))));
