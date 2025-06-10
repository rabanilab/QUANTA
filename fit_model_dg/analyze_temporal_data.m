function analyze_temporal_data()

% read parameter file
S = readstruct('param.xml');
if (isfield(S,'data_r'))
    all_data_r = cellstr(S.data_r);
else
    all_data_r = {};
end
if (isfield(S,'data_a'))
    all_data_a = cellstr(S.data_a);
else
    all_data_a = {};
end
all_data_c = cellstr(S.data_c);
mzt_class = S.mzt_class;
mzt_fit = S.mzt_fit;
kmer_range = S.kmer_range;
kmer_alpha = S.kmer_alpha;
kmer_esize = S.kmer_esize;
kmer_dir = S.kmer_dir{1};
seq_file = S.seq_file{1};
data_dir = S.data_dir{1};
class_FC = S.class_FC;
run_zfish = S.run_zfish;
run_frog = S.run_frog;
run_mouse = S.run_mouse;
run_human = S.run_human;
is_mzt = S.is_mzt;
cnt_gene_norm1 = S.cnt_gene_norm1;
cnt_gene_norm2 = S.cnt_gene_norm2;
fpkm_pref = S.fpkm_pref{1};
gene_file = S.gene_file{1};
list_file = S.list_file{1};

% define constant values
minE = 2^-3; % minimal expression flooring
minFPKM = 2; % for selecting expressed genes
minPCT = 2; % for selecting a model
minR2 = 0.5; % for selecting a good model fit
maxERR = 10; % for selecting a good model fit
class_minE = -2.5;
NORMF = 60;

nr = max(size(all_data_r));
na = max(size(all_data_a));
nc = max(size(all_data_c));

Y = importdata(gene_file);
maxN = max(nr,na);

% do not use organism specific timing
mzt_fit = [];

% ---------------------------------------------------------------------
% FPKM data
% ---------------------------------------------------------------------
load_fpkm = 0;

if (load_fpkm)
    F = importdata(list_file);

    G1 = [];
    G2 = [];
    M = [];
    P = [];
    S = [];
    T = [];
    for i = 1:size(F.textdata,1)
        fid = F.textdata{i,1};
        sid = F.textdata{i,2};
        tm = F.data(i);
        
        [data_dir '/' fid '.star.p.fpkm']
        D = importdata([data_dir '/' fid '.star.p.fpkm']);
        K = regexprep(D.textdata(2:end,1),'^gene:','');
        [~,i1,i2] = intersect(K,Y);
        G1 = K(i1);
        K = D.textdata(2:end,2);
        G2 = K(i1);
        M(:,i) = D.data(i1,1);
        P(:,i) = D.data(i1,2);
        S{i} = sid;
        T(i) = tm;
        fprintf('input: sample %s, time %.1f (n = %d)\n',sid,tm,size(i1,1));
    end
    M(M<minE) = 0;
    P(P<minE) = 0;

    u = unique(S);
    for i = 1:max(size(u))
        j = strcmp(S,u{i});
        tm = T(j);
        [~,k] = sort(tm);
        T(j) = tm(k);
        m = M(:,j);
        M(:,j) = m(:,k);
        p = P(:,j);
        P(:,j) = p(:,k);
    end
    save('data.fpkm.mat','S','T','G1','G2','M','P');
end

% ---------------------------------------------------------------------
% select datasets for gene selection and classification
% ---------------------------------------------------------------------
load('data.fpkm.mat');

j = [];
for i = 1:nc
    k = strcmp(S,all_data_c{i});
    u = unique(T(k));
    j(i) = (sum(u>=mzt_class)>0).*(sum(u<mzt_class)>0);
end
nc = sum(j);
all_data_c = all_data_c(j == 1);

fprintf('datasets selected for gene selection and classification:\n');
all_data_c

% ---------------------------------------------------------------------
% normalize expression by stable controls
% select expressed genes by maximal expression
% ---------------------------------------------------------------------
norm_fpkm = 0;

if (norm_fpkm)
    load('data.fpkm.mat');
    u = unique(S);
    MN = M;
    PN = P;

    % first step: local normalization between samples within a dataset
    % using control genes
    if (cnt_gene_norm1)
        for i = 1:max(size(u))
            j = strcmp(S,u{i});

            h = figure;
            scrsz = get(0,'ScreenSize');
            set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
            n1 = -1;
            n2 = 0;
            r = 1;
            while ((n1<n2) && (r<=15))
                x = MN(:,j) + PN(:,j);
                THR = max(prctile(x(x>minE),NORMF),2.^minFPKM);
                mnX = mean(x,2);
                stdX = std(x,1,2);
                x1 = log2(mnX);
                x2 = log2(stdX);
                k = (~isnan(x1)).*(~isnan(x2)).*(~isinf(x1)).*(~isinf(x2)) == 1;
                p = polyfit(x1(k.*(x1>2)==1),x2(k.*(x1>2)==1),1);
                c = (x1 >= log2(THR)).*(x2 <= polyval(p,x1)-1)==1;
                n1 = n2;
                n2 = sum(c);
                fprintf('norm %s, %d controls: %d\n',regexprep(u{i},'_',' '),r,n2);

                if (n2==0)
                    N = 1;
                else
                    N = mean(x(c,:));
                    N = N./mean(N);
                end
                MN(:,j) = (MN(:,j)./N);
                PN(:,j) = (PN(:,j)./N);

                subplot(3,5,r);
                if (sum(k)>0)
                    dscatter(x1(k),x2(k));
                    hold on;
                    plot(x1,polyval(p,x1),'-k','linewidth',2);
                    plot(x1,polyval(p,x1)-1,'-r','linewidth',2);
                end
                hold off;
                axis square;
                axis tight;
                xlabel('mean');
                ylabel('stdv');
                title(sprintf('norm %s, %d (n=%d, c=%d)',regexprep(u{i},'_',' '),r,sum(k),sum(c)));
                r = r+1;
            end
            saveas(h, ['results/norm.' u{i} '.jpg'],'jpg');
            close all;
        end
    end
    
    % second step: global normalization between dataset
    % by median expression
    PN2 = PN;
    MN2 = MN;
    if (cnt_gene_norm2)
        Mall = [];
        for i = 1:max(size(u))
            j = strcmp(S,u{i});
            RW = PN2(:,j) + MN2(:,j);
            Mall(1,i) = mean(RW(RW>minE));
        end

        for i = 1:max(size(u))
            j = strcmp(S,u{i});
            N = Mall(1,i)./mean(Mall(1,:));
            PN2(:,j) = PN2(:,j)./N;
            MN2(:,j) = MN2(:,j)./N;
            RW = PN2(:,j) + MN2(:,j);
            Mall(2,i) = mean(RW(RW>minE));
        end
        Mall
    end
    PN2(PN2<minE) = minE;
    MN2(MN2<minE) = minE;
    logP = log2(PN2);
    logM = log2(MN2);

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    for i = 1:max(size(u))
        subplot(2,ceil(max(size(u))/2),i);
        j = strcmp(S,u{i});
        n = sum(j);
        hold on;
        RW = M(:,j);
        RW(RW<minE) = minE;
        logRW = log2(RW);
        boxplot(logM(:,j),'widths',1,'positions',(1:n) + 0.35,'Colors','r', ...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        boxplot(logRW,'labels',T(j),'widths',1,'positions',1:n,'Colors','b',...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        hold off;
        xlabel('samples');
        ylabel('FPKM; log2');
        set(gca,'ylim',[log2(minE) 15]);
        title(regexprep(u{i},'_','-'));
    end
    saveas(h, ['results/norm.' num2str(NORMF) '.jpg'],'jpg');
        
    % select expressed genes by maximal expression
    maxE = [];
    for i = 1:nc
        k = strcmp(S,all_data_c{i});
        maxE(:,i) = max([logP(:,k) logM(:,k)],[],2);
    end
    
    k = (sum(maxE>=minFPKM,2) >= floor(nc/2)+1);
    fprintf('select by maximal expression: %d genes\n', sum(k));
    logP = logP(k,:);
    logM = logM(k,:);
    G1 = G1(k,:);
    G2 = G2(k,:);
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = log2(minE):0.25:10;
    for i = 1:nr
        subplot(2,max(maxN,3),i);
        j = strcmp(S,all_data_r{i});
        hold on;
        %violinplot(logP(:,j),T(j),'ShowData',false,'ViolinAlpha',1,'ViolinColor',[0.9 0.9 0.9],'Bandwidth',0.1)
        n = sum(j);
        boxplot(logM(:,j),'widths',1,'positions',(1:n) + 0.35,'Colors','b', ...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        boxplot(logP(:,j),'labels',T(j),'widths',1,'positions',1:n,'Colors','r',...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        hold off;
        axis tight;
        set(gca,'ylim',[log2(minE) 15]);
        ylabel('FPKM; log2');
        xlabel('samples');
        title(regexprep(all_data_r{i},'_',' '));
    end
    for i = 1:na
        subplot(2,max(maxN,3),max(maxN,3)+i);
        j = strcmp(S,all_data_a{i});
        hold on;
        %violinplot(logP(:,j),T(j),'ShowData',false,'ViolinAlpha',1,'ViolinColor',[0.9 0.9 0.9],'Bandwidth',0.1)
        n = sum(j);
        boxplot(logM(:,j),'widths',1,'positions',(1:n) + 0.35,'Colors','b', ...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        boxplot(logP(:,j),'labels',T(j),'widths',1,'positions',1:n,'Colors','r',...
            'plotstyle','compact','MedianStyle','target','Symbol','.k','OutlierSize',4,'Jitter',0.1,'whisker',2);
        hold off;
        axis tight;
        set(gca,'ylim',[-3 15]);
        ylabel('log FPKM');
        xlabel('samples');
        title(regexprep(all_data_a{i},'_',' '));
    end
    saveas(h, 'results/fpkm_hist.jpg','jpg');
    saveas(h, 'results/fpkm_hist.eps','epsc');
    
%     h = figure;
%     scrsz = get(0,'ScreenSize');
%     set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
%     for i = 1:max(size(u))
%         clf;
%         j = strcmp(S,u{i});
%         k = (T == min(T(j)));
%         vM = logM(:,j) - mean(logM(:,k),2);
%         vP = logP(:,j) - mean(logP(:,k),2);
%         s = max([vM vP],[],2) >= 1;
%         vM = vM(s,:);
%         vP = vP(s,:);
%         [~,si] = cluster_sort([vM vP]);
%         
%         subplot(1,3,1);
%         im = imagesc(vM(si,:),[-4 4]);
%         im.AlphaData = ~isnan(vM(si,:));
%         colormap(gene_colormap(1));
%         set(gca,'xtick',1:sum(j),'xticklabel',T(j),'ytick',[]);
%         xtickangle(90);
%         title(sprintf('mature mRNA (n=%d)',sum(s)));
%         subplot(1,3,2);
%         im = imagesc(vP(si,:),[-4 4]);
%         im.AlphaData = ~isnan(vP(si,:));
%         colormap(gene_colormap(1));
%         set(gca,'xtick',1:sum(j),'xticklabel',T(j),'ytick',[]);
%         xtickangle(90);
%         title(sprintf('precursor mRNA (n=%d)',sum(s)));
%         subplot(1,3,3);
%         im = imagesc(vM(si,:),[-4 4]);
%         im.AlphaData = ~isnan(vM(si,:));
%         colormap(gene_colormap(1));
%         set(gca,'xtick',1:sum(j),'xticklabel',T(j),'ytick',1:100:sum(s));
%         xtickangle(90);
%         colorbar;
%         saveas(h, ['results/fpkm_map.' u{i} '.jpg'],'jpg');
%         saveas(h, ['results/fpkm_map.' u{i} '.eps'],'epsc');
%     end

    for i = 1:max(size(u))
        j = strcmp(S,u{i});
        L = [{'ensid' 'gid'} strcat('M_', regexprep(cellstr(num2str(T(j)')),' ',''))' strcat('P_', regexprep(cellstr(num2str(T(j)')),' ',''))'];
        write_text_file(['data.normFPKM.' u{i} '.txt'],[L;[G1 G2 num2cell([logM(:,j) logP(:,j)])]]);
    end
    save('data.norm.mat','S','T','G1','G2','logM','logP');
    
    close all;
end

% ---------------------------------------------------------------------
% classification
% ---------------------------------------------------------------------
classify_genes = 0;

if (classify_genes)
    load('data.norm.mat');
    
    C = []; % 1 = M; 2 = Z; 3 = MZ
    for i = 1:nc
        k = strcmp(S,all_data_c{i});
        k0 = k.*(T < mzt_class) == 1;
        k1 = k.*(T >= mzt_class) == 1;
        
        isM = zeros(size(logM,1),1);
        isZ = zeros(size(logM,1),1);
        
        % (1) Maternal: Mature expression is high enough before MZT
        %     maxE0 = [max mature expression before MZT]
        %     T = 2x standard deviation below the mean (excluding values of -3 that skew the distribution)
        %     if maxE0 >= T, then maternal
        maxE0 = max(logM(:,k0),[],2);
        maxE0_T = mean(maxE0(maxE0>minE)) - 2*std(maxE0(maxE0>minE),1);
        isM = isM + (maxE0 >= maxE0_T) > 0;
        
        % (2) Zygotic: Mature expression after MZT is 75% or more than before
        %     maxM1 = [max mature expression after MZT]
        %        rM = [mean mature expression after MZT] / [mean mature expression before MZT]
        %     T1 = -2.5
        %     T2 = log2(0.75)
        %     if (rM >= T2) AND (maxM1 >= T1), then zygotic
        maxM1 = max(logM(:,k1),[],2);
        rM = mean(logM(:,k1),2) - mean(logM(:,k0),2);
        rM_T = log2(class_FC);
        isZ = isZ + (rM >= rM_T).*(maxM1 >= class_minE) > 0;
        
        % (3) Zygotic: Precursor expression increase by at least 25% after MZT
        %     maxP1 = [max precursor expression after MZT]
        %        rP = [max precursor expression after MZT] / [max precursor expression before MZT]
        %     T1 = -2.5
        %     T2 = log2(1.25)
        %     if (rP >= T2) AND (maxP1 >= T1), then zygotic
        if (is_mzt)
            maxP1 = max(logP(:,k1),[],2);
            rP = max(logP(:,k1),[],2) - max(logP(:,k0),[],2);
            rP_T = log2(1.25);
            isZ = isZ + (rP >= rP_T).*(maxP1 >= class_minE) > 0;
            
        % (3) Zygotic: Precursor expression after MZT is 75% or more than before
        %     maxP1 = [max mature expression after MZT]
        %        rP = [mean mature expression after MZT] / [mean mature expression before MZT]
        %     T1 = -2.5
        %     T2 = log2(0.5)
        %     if (rP >= T2) AND (maxP1 >= T1), then zygotic
        else
            maxP1 = max(logP(:,k1),[],2);
            rP = mean(logP(:,k1),2) - mean(logP(:,k0),2);
            rP_T = log2(0.5);
            isZ = isZ + (rP >= rP_T).*(maxP1 >= class_minE) > 0;
        end
        
        C(:,i) = isM + 2*isZ;
        H = cellfun(@isempty,regexp(G2,'^[Hh][1234]'));
        C(H==0,i) = 0;
        
        y = hist(C(:,i),0:3);
        fprintf('classify: %s\n', all_data_c{i});
        [{'NA' 'M' 'Z' 'MZ'}; num2cell(y); num2cell(y./sum(y))]'
    
        h = figure;
        scrsz = get(0,'ScreenSize');
        set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
        subplot(2,2,1);
        x = -3:0.2:10;
        y = hist(maxE0,x);
        hold on;
        plot(x,y./sum(y),'-k','linewidth',1.5);
        line([maxE0_T maxE0_T],[0 max(y./sum(y))],'LineStyle','-','color','r','linewidth',1);
        hold off;
        axis tight;
        xlabel('max mature expression before MZT');
        title(sprintf('n=%d above threshold (=M)',sum((maxE0 >= maxE0_T))));
        subplot(2,2,3);
        x = -2:0.1:8;
        y = hist(rP,x);
        hold on;
        plot(x,y./sum(y),'-k','linewidth',1.5);
        line([rP_T rP_T],[0 max(y./sum(y))],'LineStyle','-','color','r','linewidth',1);
        hold off;
        axis tight;
        xlabel('Precursor: after MZT / before MZT; log2');
        title(sprintf('n=%d above threshold (=Z)',sum((rP >= rP_T))));
        subplot(2,2,4);
        x = -8:0.5:8;
        y = hist(rM,x);
        hold on;
        plot(x,y./sum(y),'-k','linewidth',1.5);
        line([rM_T rM_T],[0 max(y./sum(y))],'LineStyle','-','color','r','linewidth',1);
        hold off;
        axis tight;
        xlabel('Mature: after MZT / before MZT; log2');
        title(sprintf('n=%d above threshold (=Z)',sum((rP >= rP_T))));
        saveas(h, ['results/class.' all_data_c{i} '.jpg'],'jpg');       
    end
    
    Cx = cell(size(C));
    Cx(C==1) = {'M'};
    Cx(C==2) = {'Z'};
    Cx(C==3) = {'MZ'};
    Ce = repmat({'NA'},size(Cx,1),1);
    for i = 1:nc
        k = strcmp(S,all_data_c{i});
        titles = ['class'; 'ensid'; 'gid'; regexprep(cellstr(strcat('P_',num2str(T(k)'))),' ',''); regexprep(cellstr(strcat('M_',num2str(T(k)'))),' ','')]';
        write_text_file([fpkm_pref all_data_c{i} '.txt'], [titles;[Cx(:,i) G1 G2 num2cell([logP(:,k) logM(:,k)])]]);
    end
    u = setdiff(unique(S),all_data_c);
    for i = 1:size(u,2)
        k = strcmp(S,u{i});
        titles = ['class'; 'ensid'; 'gid'; regexprep(cellstr(strcat('P_',num2str(T(k)'))),' ',''); regexprep(cellstr(strcat('M_',num2str(T(k)'))),' ','')]';
        write_text_file([fpkm_pref u{i} '.txt'], [titles;[Ce G1 G2 num2cell([logP(:,k) logM(:,k)])]]);
    end
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    elim = [-3 16];
    for j = 1:nc
        subplot(2,max(ceil(nc/2),3),j);
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
        k0 = ismember(tid,Y);
        dscatter(max(D(k0,:),[],2),mean(D(k0,iM),2),'MSIZE',20);
        xlabel('max');
        ylabel('mean');
        set(gca,'xlim',elim,'ylim',elim);
        axis square;
        title(sprintf('%s (n=%d, max:%d, mean:%d)',all_data_c{j},sum(k0),sum(max(D(k0,:),[],2)>=2),sum(mean(D(k0,iM),2)>0)));
        hold on;
        line(elim,elim,'LineStyle','-','color','k','linewidth',1);
        line([2 2],elim,'LineStyle','-','color','k','linewidth',1);
        line(elim,[0 0],'LineStyle','-','color','k','linewidth',1);
        hold off;
    end
    saveas(h, 'results/expression.jpg','jpg');

    save('data.class.mat','S','T','G1','G2','logM','logP','C');
    close all;
end

% ---------------------------------------------------------------------
% combined classification
% ---------------------------------------------------------------------
class_plot = 0;

if (class_plot == 1)
    S = {'M' 'Z' 'MZ'};
    allid = [];
    for j = 1:nc
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
        k0 = ismember(tid,Y);
        k0 = k0.*(max(D(:,iM),[],2)>=minFPKM)==1;
        sum(k0)
        allid = [allid;tid(k0)];

        [u,~,t] = unique(gc(k0));
        n = accumarray(t,1);
        all_data_c{j}
        [u num2cell([n n./sum(n)])]
    end
    [u,~,t] = unique(allid);
    c0 = accumarray(t,1);
    i0 = (c0 >= floor(nc/2)+1);
    allid = u(i0);
    fprintf('gene classification: %d genes\n', max(size(allid)));
    
    mtid = [];
    ztid = [];
    for j = 1:nc
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
        k0 = ismember(tid,allid);
        
        k = strcmp(gc,'M').*k0==1;
        mtid = [mtid;tid(k)];
        k = strcmp(gc,'Z').*k0==1;
        ztid = [ztid;tid(k)];
    end
    [u,~,t] = unique(mtid);
    cm = accumarray(t,1);
    im = (cm >= floor(nc/2)+1);
    mtid = u(im);
    
    [u,~,t] = unique(ztid);
    cz = accumarray(t,1);
    iz = (cz >= floor(nc/2)+1);
    ztid = u(iz);
    
    xc = 3*ones(size(tid));
    xc(ismember(tid,mtid)) = 1;
    xc(ismember(tid,ztid)) = 2;
    
    fprintf('gene classification:\n');
    y = hist(xc,1:3);
    [S' num2cell([y' y'./sum(y)])]
    save('gene_classification.mat','mtid','ztid','allid');
    L = {'ensid' 'gid' 'class'};
    write_text_file('gene_classification.txt',[L;[tid gid num2cell(xc)]]);
    
    h = figure;
    b = bar(y./sum(y));
    axis tight;
    ylabel('fraction of genes');
    set(gca,'xticklabel',S,'ylim',[0 1],'fontsize',14);
    for i = 1:3
        text(i-0.3,y(i)/sum(y)-0.05,sprintf('n=%d (%.1f%%)',y(i),100*y(i)/sum(y)),'fontsize',14);
    end
    b.FaceColor = 'flat';
    b.CData(1,:) = [1 0 0];
    b.CData(2,:) = [0 0 1];
    b.CData(3,:) = [0.8 0.8 0.8];
    saveas(h, 'results/class.jpg','jpg');
    saveas(h, 'results/class.eps','epsc');
    close all;

    if (nc > 1)
        mid = [];
        zid = [];
        mzid = [];
        for j = 1:nc
            [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
            k0 = ismember(tid,allid);
            
            k = strcmp(gc,'M').*k0==1;
            mid = [mid;tid(k)];
            k = strcmp(gc,'Z').*k0==1;
            zid = [zid;tid(k)];
            k = strcmp(gc,'MZ').*k0==1;
            mzid = [mzid;tid(k)];
        end
        [u,~,t] = unique(mid);
        cm = accumarray(t,1);
        ym = hist(cm,1:nc);
        
        [u,~,t] = unique(zid);
        cz = accumarray(t,1);
        yz = hist(cz,1:nc);
        
        [u,~,t] = unique(mzid);
        cmz = accumarray(t,1);
        ymz = hist(cmz,1:nc);
        
        fprintf('gene classification: compare between datasets\n');
        C = [ym;yz;ymz];
        C = C(:,2:nc);
        T = {' ' 'all'};
        for i = 2:nc
            T = [T sprintf('%d/%d',i,nc)];
        end
        [[T;[S' num2cell([y' C])]]; ...
            [{'all';'all(%)'} num2cell([sum([y' C]);sum([y' C])./sum(y)])]]
        
        h = figure;
        f = sum(C);
        f = [sum(y)-sum(f) f];
        b = bar(f./sum(f));
        axis tight;
        ylabel('fraction of genes');
        T = {'rest'};
        for i = 2:nc
            T = [T sprintf('%d/%d',i,nc)];
        end
        set(gca,'xticklabel',T,'ylim',[0 1],'fontsize',14);
        for i = 1:nc
            text(i-0.3,f(i)./sum(f)-0.05,sprintf('n=%d (%.1f%%)',f(i),100*f(i)./sum(f)),'fontsize',14);
        end
        b.FaceColor = 'flat';
        for i = 1:nc
            b.CData(i,:) = [0.8 0.8 0.8];
        end
        saveas(h, 'results/class.compare.jpg','jpg');
        saveas(h, 'results/class.compare.eps','epsc');
        close all;
    end
    
    if (run_zfish)
        f = fopen('/Users/mrabani/GoogleDrive/Rabani_Lab/Projects/project_MZT_transcriptomics/0_zebrafish_Nov2022/1-s2.0-S2211124723000815-mmc4.txt');
        X1 = textscan(f,'%s %s');
        fclose(f);
        f = fopen('/Users/mrabani/GoogleDrive/Rabani_Lab/Projects/project_MZT_transcriptomics/0_zebrafish_Nov2022/model_classification_scSLAM.txt');
        X2 = textscan(f,'%s %s');
        fclose(f);
       
        [~,i1,i2] = intersect(X1{1},X2{1});
        yg = X1{1}(i1);
        class = [X1{2}(i1) X2{2}(i2(i2>0))];
        yc = zeros(size(class));
        yc(strcmp(class,'M')) = 1;
        yc(strcmp(class,'Z')) = 2;
        yc(strcmp(class,'MZ')) = 3;
        
        xg = gid;
        [~,i1,i2] = intersect(xg,yg);
        zg = xg(i1); 
        zc = [xc(i1,:) yc(i2(i2>0),:)];
        
        % [this study] [Bhat et. al.] [Fishman et. al]
        [u,~,t] = unique(zc,'rows');
        n = accumarray(t,1);
        num2cell(sortrows([u n n./sum(n)]))
        
        h = figure;
        f = [];
        f(1) = sum(n((u(:,1)==u(:,2)).*(u(:,2)==u(:,3))==1));
        f(2) = sum(n((u(:,1)==u(:,2)).*(u(:,1)~=u(:,3))==1));
        f(3) = sum(n((u(:,1)~=u(:,2)).*(u(:,1)==u(:,3))==1));
        f(4) = sum(n((u(:,1)~=u(:,2)).*(u(:,2)==u(:,3))==1));
        f = [f sum(n)-sum(f)];
        b = bar(f./sum(f));
        axis tight;
        ylabel('fraction of genes');
        set(gca,'xticklabel',{'3/3' '2/3 (bhat)' '2/3 (fishman)' '2/3 (other)' 'rest'},'ylim',[0 1],'fontsize',14);
        b.FaceColor = 'flat';
        for i = 1:5
            text(i-0.3,f(i)./sum(f)-0.05,sprintf('n=%d (%.1f%%)',f(i),100*f(i)./sum(f)),'fontsize',14);
            b.CData(i,:) = [0.8 0.8 0.8];
        end
        saveas(h, 'results/class.compare.zf.jpg','jpg');
        saveas(h, 'results/class.compare.zf.eps','epsc');
        close all;
        
        for i = 1:3
            k = zc(:,1)==i;
            [u,~,t] = unique(zc(k,:),'rows');
            n = accumarray(t,1);
            num2cell(sortrows([u n n./sum(n)],4))
        end
    end
end

% ---------------------------------------------------------------------
% maternal genes
% ---------------------------------------------------------------------
load('gene_classification.mat','mtid');
nmtid = size(mtid,1);

% ---------------------------------------------------------------------
% expression data: stable genes
% ---------------------------------------------------------------------
plot_controls = 0;

if (plot_controls)
    load('gene_classification.mat','allid');
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    mD = cell(1,nc);
    sD = cell(1,nc);
    for j = 1:nc
        all_data_c{j}
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
        [k1,k2] = ismember(allid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        stdM = std(y,1,2);
        mnM = mean(y,2);
        mD{j} = nan(size(allid,1),1);
        mD{j}(k1,:) = mnM;
        sD{j} = nan(size(allid,1),1);
        sD{j}(k1,:) = stdM;
        
        subplot(2,max(ceil(nc/2),3),j);
        dscatter(mnM,stdM);
        hold on;
        line([0 15],[1 1],'LineStyle','-','color','k','linewidth',1);
        line([9 9],[0 12],'LineStyle','-','color','k','linewidth',1);
        hold off;
        axis square;
        set(gca,'ylim',[0 8],'xlim',[0 15]);
        xlabel('mean');
        ylabel('stdv');
        k = (mnM>=9).*(stdM<=1)==1;
        title(sprintf('%s (%d)',all_data_c{j},sum(k)));
    end
    saveas(h, 'results/controls.jpg','jpg');
    close all;
            
    mkdir('results/expr_c');
    mW = [mD{:}];
    sW = [sD{:}];
    k = find((sum(mW >= 5,2) == nc).*(sum(sW <= 2,2) == nc)==1);
    s = sortrows([mean(sW(k,:),2) mean(mW(k,:),2) k]);
    last = min(size(s,1),10);
    plot_ids(allid(s(1:last,3)),all_data_a,all_data_r,'results/expr_c/',0,0,mzt_fit,0,fpkm_pref);
end

% ---------------------------------------------------------------------
% expression data: plot examples
% ---------------------------------------------------------------------
plot_examples = 0;

if (plot_examples)
    load('gene_classification.mat','mtid','ztid','allid');
    mztid = setdiff(allid,mtid);
    mztid = setdiff(mztid,ztid);

    mD = cell(1,nc);
    for j = 1:nc
        all_data_c{j}
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_c{j} '.txt']);
        [k1,k2] = ismember(allid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        mnM = mean(y,2);
        mD{j} = nan(size(allid,1),1);
        mD{j}(k1,:) = mnM;
    end
    mW = [mD{:}];
    
    mkdir('results/expr_select');
    k = find(ismember(allid,mtid));
    if (~isempty(k))
        s = sortrows([mean(mW(k,:),2) k],-1);
        last = min(size(s,1),15);
        plot_ids(allid(s(1:last,2)),all_data_a,all_data_r,'results/expr_select/mx_',1,0,mzt_fit,0,fpkm_pref);
    end
    k = find(ismember(allid,ztid));
    if (~isempty(k))
        s = sortrows([mean(mW(k,:),2) k],-1);
        last = min(size(s,1),5);
        plot_ids(allid(s(1:last,2)),all_data_a,all_data_r,'results/expr_select/zx_',0,0,mzt_fit,0,fpkm_pref);
    end
    k = find(ismember(allid,mztid));
    if (~isempty(k))
        s = sortrows([mean(mW(k,:),2) k],-1);
        last = min(size(s,1),5);
        plot_ids(allid(s(1:last,2)),all_data_a,all_data_r,'results/expr_select/mzx_',0,0,mzt_fit,0,fpkm_pref);
    end
    
    if (run_zfish)
        all_ids = {'dazl' 'buc' 'slbp2' 'btg4' 'h1m' 'wee2' 'thy1' 'tdrd5' 'tom1' 'lrmp' ...
            'cdc34a' 'kpna7' 'ca15b' 'cntd2' 'birc5b' 'org' 'aktip' 'gdf3' 'cpeb1b' 'ddx4' ...
            'ccna1' 'slc16a3' 'zp3.2' 'nsdhl' 'pja2' 'impdh2' 'lman1' 'maza' 'avd' ...
            'pttg1' 'pdlim2' 'mctp2b' 'mcm3l' 'mgat2' 'washc2c' 'fbxl14a' 'tspan36' 'necap1' ...
            'si:dkeyp-114g9.1' 'zgc:171776' 'zgc:110239' 'zgc:162879'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/m_',1,1,mzt_fit,0,fpkm_pref);
        
        all_ids = {'tbxta' 'eve1' 'cdx4' 'ved' 'gsc' 'sox3' 'nog1' 'sox32' 'spaca4l' 'krt4' 'krt8'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/z_',0,1,mzt_fit,0,fpkm_pref);
        
        all_ids = {'scube2' 'spry2' 'pou5f3' 'sox11b' 'nanog' 'pfn1' 'cth1' 'xbp1' 'fscn1a' 'ccng1'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/mz_',0,1,mzt_fit,0,fpkm_pref);
        
        all_ids = {'smad7' 'zic2b' 'galnt2' 'foxk1' 'trim25l' 'me3' 'neo1b' 'dynlrb1'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/u_',0,1,mzt_fit,0,fpkm_pref);
        
        % eps plots
        all_ids = {'dazl' 'buc' 'lhx8a' 'slbp2' 'btg4' 'wee2'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/x_',1,1,mzt_fit,1,fpkm_pref);
        
        all_ids = {'eve1' 'fgfr4' 'pum3' 'xbp1' 'rps26' 'tubb4b'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/x_',0,1,mzt_fit,1,fpkm_pref);
    end
    if (run_frog)
        all_ids = {'velo1' 'h1-8' 'dazl' 'btg4' 'wee2' 'tdrd5' 'tom1' 'lrmp' ...
            'kpna7' 'cntd2' 'aktip' 'ccna1' 'zp3.2' 'nsdhl' 'necap1'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/m_',1,1,mzt_fit,0,fpkm_pref);
    end        
    if (run_mouse)
        all_ids = {'Btg4' 'Wee2' 'Dazl' 'Slbp' 'Zar1' 'Dazap2' 'Ddx11' ...
            'Ddx28' 'Kif3a' 'Axin2' 'Dcp2' 'Dicer1' 'Galm'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/m_',1,1,mzt_fit,0,fpkm_pref);
    end
    if (run_human)
        all_ids = {'BTG4' 'WEE2' 'DAZL' 'TCL1A' 'H1-8' 'TUBB8' 'OOSP2' ...
            'BCAR4' 'BOD1' 'ACTL8' 'OTX2' 'MYOCOS'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/m_',1,1,mzt_fit,0,fpkm_pref);
        all_ids = {'EEF1A1' 'GAPDH' 'SLBP' 'CDC6' 'CDK1' 'ACE2' 'HHEX'};
        plot_ids(all_ids,all_data_a,all_data_r,'results/expr_select/x_',0,1,mzt_fit,0,fpkm_pref);
    end
end

% ---------------------------------------------------------------------
% model fits
% ---------------------------------------------------------------------
model_fit = 0;

if (model_fit == 1)
    load('gene_classification.mat','mtid');

    % fit 1-rate model
    if (isempty(mzt_fit))
        run_fit(mtid,[all_data_r all_data_a],[],fpkm_pref);
    else
        run_fit(mtid,[all_data_r all_data_a],mzt_fit(2:3),fpkm_pref);
    end

    % fit 2-rate model
    if (isempty(mzt_fit))
        run_fit_2p(mtid,[all_data_r all_data_a],[],fpkm_pref);
    else
        run_fit_2p(mtid,[all_data_r all_data_a],mzt_fit,fpkm_pref);
    end
end

% ---------------------------------------------------------------------
% STD
% ---------------------------------------------------------------------
% s = [];
% for j = 1:nr
%     [D,~,iM,~,~,tid] = load_data([fpkm_pref all_data_r{j} '.txt']);
%     [~,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
%     y = D(k2(k2>0),iM);
%     y = y - y(:,1);
%     s(j) = sum(sum((y(:,1:end-1)-y(:,2:end)).^2))/(size(y,1)*(size(y,2)-1));
% end
% for j = 1:na
%     [D,~,iM,~,~,tid] = load_data([fpkm_pref all_data_a{j} '.txt']);
%     [~,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
%     y = D(k2(k2>0),iM);
%     y = y - y(:,1);
%     s(nr+j) = sum(sum((y(:,1:end-1)-y(:,2:end)).^2))/(size(y,1)*(size(y,2)-1));
% end
% [[all_data_r'; all_data_a'] num2cell(s')]
% mean(s)

% ---------------------------------------------------------------------
% test model fit
% ---------------------------------------------------------------------
test_model_fit = 0;

if (test_model_fit)
    mkdir('results/model_fit');
    
    % 1-rate model results
    mP1 = cell(1,nr+na); % P = <logX0,dg,t0,t1>
    mR1 = cell(1,nr+na);
    mE1 = cell(1,nr+na);
    for j = 1:nr
        load(['model_param.1p.' all_data_r{j} '.mat'],'mtid','mP','mR','mE');
        mP1{j} = mP;
        mR1{j} = mR;
        mE1{j} = mE;
    end
    for j = 1:na
        load(['model_param.1p.' all_data_a{j} '.mat'],'mtid','mP','mR','mE');
        mP1{nr+j} = mP;
        mR1{nr+j} = mR;
        mE1{nr+j} = mE;
    end
    h = fit_model_plot_param(mP1,[all_data_r all_data_a]);
    saveas(h, 'results/model_fit/test.1p.param.jpg','jpg');
    %saveas(h, 'results/model_fit/test.1p.param.svg','svg');
    model_test_fit([all_data_r all_data_a],mtid,mP1,mE1,mR1,1,'results/model_fit/test.1p',minR2,maxERR,[],fpkm_pref);
    close all;

    % 2-rate model results
    mP2 = cell(1,nr+na); % P = <logX0,da,dg,t0,t1>
    mR2 = cell(1,nr+na);
    mE2 = cell(1,nr+na);
    for j = 1:nr
        load(['model_param.2p.' all_data_r{j} '.mat'],'mtid','mP','mR','mE');
        mP2{j} = mP;
        mR2{j} = mR;
        mE2{j} = mE;
    end
    for j = 1:na
        load(['model_param.2p.' all_data_a{j} '.mat'],'mtid','mP','mR','mE');
        mP2{nr+j} = mP;
        mR2{nr+j} = mR;
        mE2{nr+j} = mE;
    end
    h = fit_model_2p_plot_param(mP2,[all_data_r all_data_a]);
    saveas(h, 'results/model_fit/test.2p.param.jpg','jpg');
    model_test_fit([all_data_r all_data_a],mtid,mP2,mE2,mR2,2,'results/model_fit/test.2p',minR2,maxERR,[],fpkm_pref);

    % number of geness with R2 above threshold by each model
    xr = [];
    if (nr > 0)
        for j = 1:nr
            xr(j,1) = sum(mR1{j}>minR2);
            xr(j,2) = sum(mR2{j}>minR2);
        end
        xr = [xr;sum(xr,1)]./[repmat(nmtid,nr,1);nmtid*nr];
        [[all_data_r 'all']' num2cell(100*[xr xr(:,2)-xr(:,1)])]
    end
    xa = [];
    if (na>0)
        for j = nr+1:nr+na
            xa(j-nr,1) = sum(mR1{j}>minR2);
            xa(j-nr,2) = sum(mR2{j}>minR2);
        end
        xa = [xa;sum(xa,1)]./[repmat(nmtid,na,1);nmtid*na];
        [[all_data_a 'all']' num2cell(100*[xa xa(:,2)-xa(:,1)])]
    end

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    if (nr > 0)
        subplot(1,2,1);
        bar(xr);
        for j = 1:size(xr,1)
            text(j+0.1,xr(j,1),sprintf('%.0f%%',100*xr(j,1)),'FontSize',18);
            text(j-0.1,xr(j,2),sprintf('%.0f%%',100*xr(j,2)),'FontSize',18);
        end
        set(gca,'xticklabel',[all_data_r 'all'],'ylim',[0 1],'fontsize',16);
        title('total RNA');
        ylabel(sprintf('r2 > %.1f',minR2));
        legend({'degradation model' 'polyA model'},'location','northoutside','box','off');
    end
    if (na > 0)
        subplot(1,2,2);
        bar(xa);
        for j = 1:size(xa,1)
            text(j+0.1,xa(j,1),sprintf('%.0f%%',100*xa(j,1)),'FontSize',18);
            text(j-0.1,xa(j,2),sprintf('%.0f%%',100*xa(j,2)),'FontSize',18);
        end
        set(gca,'xticklabel',[all_data_a 'all'],'ylim',[0 1],'fontsize',16);
        title('polyA+ RNA');
        ylabel(sprintf('r2 > %.1f',minR2));
        legend({'degradation model' 'polyA model'},'location','northoutside','box','off');
    end
    saveas(h,'results/model_fit/select.r2.jpg','jpg');
    saveas(h,'results/model_fit/select.r2.eps','epsc');
        
    % likelihood ratio test
    h = model_llr(all_data_r, all_data_a, mtid, fpkm_pref);
    saveas(h, 'results/model_fit/select.llr.jpg','jpg');
    saveas(h, 'results/model_fit/select.llr.eps','epsc');    

    % goodness of fit test
    h = model_gof(all_data_r, all_data_a, mtid, fpkm_pref);
    saveas(h, 'results/model_fit/select.gof.jpg','jpg');
    saveas(h, 'results/model_fit/select.gof.eps','epsc');
    
    close all;

    % ribo-depleted: select a model per dataset
    Pr = [];
    Er = [];
    Rr = [];
    Jr = zeros(1,nr);
    for j = 1:nr
        % if (xr(j,2)-xr(j,1)>(minPCT/100)) % 2-rate model
        %     Pr{j} = mP2{j};
        %     Er{j} = mE2{j};
        %     Rr{j} = mR2{j};
        %     fprintf('dataset %s (total-RNA): select 2-rate model\n',all_data_r{j});
        %     Jr(j) = 1;
        % else % 1-rate model
            Pr{j} = [mP1{j}(:,1) zeros(size(mP1{j},1),1) mP1{j}(:,2:4)];
            Er{j} = mE1{j};
            Rr{j} = mR1{j};
        % end
    end

    % polyA+: select a model per dataset
    Pa = [];
    Ea = [];
    Ra = [];
    Ja = zeros(1,na);
    for j = 1:na
        if (xa(j,2)-xa(j,1)>(minPCT/100)) % 2-rate model
            Pa{j} = mP2{nr+j};
            Ea{j} = mE2{nr+j};
            Ra{j} = mR2{nr+j};
            fprintf('dataset %s (polyA-RNA): select 2-rate model\n',all_data_a{j});
            Ja(j) = 1;
        else % 1-rate model
            Pa{j} = [mP1{nr+j}(:,1) zeros(size(mP1{nr+j},1),1) mP1{nr+j}(:,2:4)];
            Ea{j} = mE1{nr+j};
            Ra{j} = mR1{nr+j};
        end
    end
    save('maternal_param.mat','mtid','Pa','Pr','Ra','Rr','Ea','Er','Jr','Ja');

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    if (nr > 0)
        subplot(1,2,1);
        bar(Jr);
        set(gca,'xticklabel',regexprep(all_data_r,'_',' '),'ylim',[0 1],'fontsize',16);
        title('total RNA');
        ylabel(sprintf('selected 2-rate model (>%.1f%% of genes improved)',minPCT));
    end
    if (na > 0)
        subplot(1,2,2);
        bar(Ja);
        set(gca,'xticklabel',regexprep(all_data_a,'_',' '),'ylim',[0 1],'fontsize',16);
        title('polyA+ RNA');
        ylabel(sprintf('selected 2-rate model (>%.1f%% of genes improved)',minPCT));
    end
    saveas(h,'results/model_fit/select.jpg','jpg');
    saveas(h,'results/model_fit/select.eps','epsc');

    h = fit_model_2p_plot_param([Pr Pa],[all_data_r all_data_a]);
    saveas(h, 'results/model_fit/test.final.param.jpg','jpg');
    if ((na>0)&&(nr>0))
        h = compare_pa_rz(mtid,all_data_r,Pr,all_data_a,Pa);
        saveas(h, 'results/model_fit/compare.a.r.jpg','jpg');
    end

    close all;
end

% ---------------------------------------------------------------------
% average model fits across datasets
% ---------------------------------------------------------------------
model_average = 0;

if (model_average == 1)
    load('maternal_param.mat','mtid','Pa','Pr','Ra','Rr','Ea','Er','Ja','Jr');

    % correlations
    E2r = [];
    R2r = [];
    for j = 1:nr
        R2r(:,j) = Rr{j};
        E2r(:,j) = Er{j};
    end
    E2a = [];
    R2a = [];
    for j = 1:na
        R2a(:,j) = Ra{j};
        E2a(:,j) = Ea{j};
    end
    Ir = (R2r>=minR2) + (E2r<=maxERR) > 0;
    Ia = (R2a>=minR2) + (E2a<=maxERR) > 0;
    s0 = sum([Ir Ia]==0);
    s1 = sum([Ir Ia]==1);
    [[all_data_r all_data_a]' num2cell([s1' s0' (s0./(s0+s1))'])]
    
    % degradation rate
    Mr_dg = nan(nmtid,1);
    for j = 1:nr
        Mr_dg(:,j) = Pr{j}(:,3);
    end
    Mr_dg(Ir==0) = NaN;
    Ma_dg = nan(nmtid,1);
    for j = 1:na
        Ma_dg(:,j) = Pa{j}(:,3);
    end
    Ma_dg(Ia==0) = NaN;
    fprintf('degradation rates:\n');
    corr([Mr_dg Ma_dg],'rows','pairwise')
    R(:,1) = exp(nanmean(log([Mr_dg Ma_dg]),2));
    
    % half-life
    Mr_hl = log(2)./Mr_dg;
    Ma_hl = log(2)./Ma_dg;
    fprintf('half-lives:\n');
    corr([Mr_hl Ma_hl],'rows','pairwise')
    R(:,2) = exp(nanmean(log([Mr_hl Ma_hl]),2));

    % onset time
    Mr_t0 = nan(nmtid,1);
    for j = 1:nr
        Mr_t0(:,j) = Pr{j}(:,4);
    end
    Mr_t0(Ir==0) = NaN;
    Ma_t0 = nan(nmtid,1);
    for j = 1:na
        Ma_t0(:,j) = Pa{j}(:,4);
    end
    Ma_t0(Ia==0) = NaN;
    fprintf('onset times:\n');
    corr([Mr_t0 Ma_t0],'rows','pairwise')
    R(:,3) = nanmean([Mr_t0 Ma_t0],2);

    % dA rate
    Ma_da = nan(nmtid,1);
    if (sum(Ja)>0)
        for j = 1:na
            Ma_da(:,j) = Pa{j}(:,2);
        end
        Ma_da(Ia==0) = NaN;
        Ma_da(:,Ja==0) = NaN; % 2-rate model was not selected
        fprintf('deadenylation rates:\n');
        corr(Ma_da,'rows','pairwise')
        R(:,4) = nanmean(Ma_da,2);
    else
        R(:,4) = nan(nmtid,1);
    end

    % X0 levels
    Mr_x0 = nan(nmtid,1);
    for j = 1:nr
        Mr_x0(:,j) = Pr{j}(:,1);
    end
    Mr_x0(Ir==0) = NaN;
    Ma_x0 = nan(nmtid,1);
    for j = 1:na
        Ma_x0(:,j) = Pa{j}(:,1);
    end
    Ma_x0(Ia==0) = NaN;
    fprintf('initial expression:\n');
    corr([Mr_x0 Ma_x0],'rows','pairwise')
    R(:,5) = nanmean([Mr_x0 Ma_x0],2);

    % Temporal expression pA/RT ratio
    if ((nr>0)*(na>0)==1)
        mxT = 0;
        for j = 1:nr
            [~,t] = load_data([fpkm_pref all_data_r{j} '.txt']);
            mxT = max(max(t),mxT);
        end
        for j = 1:na
            [~,t] = load_data([fpkm_pref all_data_a{j} '.txt']);
            mxT = max(max(t),mxT);
        end
        xT = 0:mxT;
        nTx = size(xT,2);
        if (nTx > 15)
            d = round(nTx/15);
            xT = xT(:,1:d:nTx);
        end

        Xr = zeros(nmtid,max(size(xT))); % ribo RNA
        Nr = zeros(nmtid,1);
        for j = 1:nr
            X = [];
            for i = 1:nmtid
                X(i,:) = dg_eval_model_2p(xT,Pr{j}(i,:));
            end
            Xr = Xr + X;
            Nr = Nr + (~isnan(Pr{j}(:,1)));
        end
        Xr = Xr./Nr;
        Xr(isnan(Xr)) = -10;
        
        Xa = zeros(nmtid,max(size(xT))); % polyA RNA
        Na = zeros(nmtid,1);
        for j = 1:na
            X = [];
            for i = 1:nmtid
                X(i,:) = dg_eval_model_2p(xT,Pa{j}(i,:));
            end
            Xa = Xa + X;
            Na = Na + (~isnan(Pa{j}(:,1)));
        end
        Xa = Xa./Na;
        Xa(isnan(Xa)) = -10;

        Rx = Xa-Xr;
        Rx_id = regexprep(cellstr(strcat(num2str(xT'),'hr')),' ','')';
        corr(Rx)

        h = plot_ATR(xT,Xa,Xr,Rx);
        saveas(h,'results/param.ratio.jpg','jpg');
    else
        Xa = [];
        Xr = [];
        Rx = [];
        Rx_id = [];
    end

    save('maternal_param.avg.mat','mtid','Ir','Ia','Mr_dg','Ma_dg',...
        'Mr_hl','Ma_hl','Ma_da','Mr_t0','Ma_t0','Mr_x0','Ma_x0',...
        'R','Rx','Rx_id','Xa','Xr');
    [~,~,~,~,gid,tid] = load_data([fpkm_pref all_data_c{1} '.txt']);
    [~,i1,i2] = intersect(mtid,tid); % [mtid(i1) tid(i2(i2>0))]
    L = {'ensid' 'gid' 'dg' 'hl' 't0' 'da' 'x0'};
    write_text_file('maternal_param.avg.txt', [L;[mtid(i1,:) gid(i2(i2>0)) num2cell(R(i1,:))]]);
    L = [{'ensid' 'gid'} all_data_r all_data_a];
    if (nr==0)
        write_text_file('maternal_param.HL.txt', [L;[mtid(i1,:) gid(i2(i2>0)) num2cell(Ma_hl(i1,:))]]);
    elseif (na==0)
        write_text_file('maternal_param.HL.txt', [L;[mtid(i1,:) gid(i2(i2>0)) num2cell(Mr_hl(i1,:))]]);
    else
        write_text_file('maternal_param.HL.txt', [L;[mtid(i1,:) gid(i2(i2>0)) num2cell([Mr_hl(i1,:) Ma_hl(i1,:)])]]);
    end
    
    % plot average parameter distributions
    X = cell(1,1);
    X{1} = [R(:,5) R(:,4) R(:,1) R(:,3) nan(nmtid,1)];
    h = fit_model_2p_plot_param(X,{'avg'},{'id' 'x0' 'da' 'dg' 't0' 'n/a'},[-3 -3 0.1 0 -1],[12 3 3 11 1]);
    xlabel('n/a');
    saveas(h,'results/param.avg.jpg','jpg');
    saveas(h,'results/param.avg.eps','epsc');

    load('../1_zebrafish_Nov2023/maternal_param.avg.mat','R');
    X{2} = [R(:,5) R(:,4) R(:,1) R(:,3) nan(size(R,1),1)];
    load('maternal_param.avg.mat','R');
    h = fit_model_2p_plot_param(X,{'avg' 'zfish'},{'id' 'x0' 'da' 'dg' 't0' 'n/a'},[-3 -3 0.1 0 -1],[12 3 3 11 1]);
    xlabel('n/a');
    saveas(h,'results/param.avg.zfish.jpg','jpg');
    saveas(h,'results/param.avg.zfish.eps','epsc');
    
    %load('../../project_utrseq/utrseq_dT/Results/analyze_20180129/maternal_param.avg.mat','R0','R40');
    %X{3} = [R0(:,5) R0(:,4) R0(:,1) R0(:,3) nan(size(R0,1),1)];
    %X{4} = [R40(:,5) R40(:,4) R40(:,1) R40(:,3) nan(size(R40,1),1)];
    %h = fit_model_2p_plot_param(X,{'avg' 'zfish' 'A-' 'A+'},{'id' 'x0' 'da' 'dg' 't0' 'n/a'},[-3 -3 0.1 0 -1],[12 3 3 11 1]);
    %xlabel('n/a');
    %saveas(h,'results/param.avg.all.jpg','jpg');
    %saveas(h,'results/param.avg.all.eps','epsc');

    Y = cell(1,2);
    Y{1} = [nanmean(Mr_x0,2) nan(nmtid,1) nanmean(Mr_dg,2) nanmean(Mr_t0,2) nan(nmtid,1)];
    Y{2} = [nanmean(Ma_x0,2) nanmean(Ma_da,2) nanmean(Ma_dg,2) nanmean(Ma_t0,2) nan(nmtid,1)];
    h = fit_model_2p_plot_param(Y,{'total' 'polyA'},{'id' 'x0' 'da' 'dg' 't0' 'n/a'},[-3 -3 0.1 0 -1],[12 3 3 11 1]);
    xlabel('n/a');
    saveas(h,'results/param.avg.ar.jpg','jpg');
    saveas(h,'results/param.avg.ar.eps','epsc');

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    subplot(1,3,1);
    fprintf('all (n=%d):\n',nmtid); 
    S = [];
    for i = [1 3 4]
        S = [S;[nanmean([Y{1}(:,i) Y{2}(:,i)]) 2.^(nanmean(Y{1}(:,i))-nanmean(Y{2}(:,i)))]];
    end
    [{'id' 'ribo' 'polya' 'fold-diff'};[{'x0' 'dg' 't0'}' num2cell(round(S,2))]]
    bar(S(:,1:2));
    set(gca,'ylim',[0 6],'xticklabel',{'x0' 'dg' 't0'},'fontsize',16);
    ylabel('mean value');
    title(sprintf('all (n=%d):\n',nmtid));
    legend({'ribo' 'polya'},'location','northoutside','box','off');
    subplot(1,3,2);
    k = Y{2}(:,2) < -0.15;%-0.2;
    fprintf('deA < 0, polyadenylation (n=%d):\n',sum(k));
    S = [];
    for i = [1 3 4]
        S = [S;[nanmean([Y{1}(k,i) Y{2}(k,i)]) 2.^(nanmean(Y{1}(k,i))-nanmean(Y{2}(k,i)))]];
    end
    [{'id' 'ribo' 'polya' 'fold-diff'};[{'x0' 'dg' 't0'}' num2cell(round(S,2))]]
    bar(S(:,1:2));
    set(gca,'ylim',[0 6],'xticklabel',{'x0' 'dg' 't0'},'fontsize',16);
    ylabel('mean value');
    title(sprintf('deA < 0, polyadenylation (n=%d):\n',sum(k)));
    legend({'ribo' 'polya'},'location','northoutside','box','off');
    subplot(1,3,3);
    k = Y{2}(:,2) > 0.15;%0.2;
    fprintf('deA > 0, deadenylation (n=%d):\n',sum(k));
    S = [];
    for i = [1 3 4]
        S = [S;[nanmean([Y{1}(k,i) Y{2}(k,i)]) 2.^(nanmean(Y{1}(k,i))-nanmean(Y{2}(k,i)))]];
    end
    [{'id' 'ribo' 'polya' 'fold-diff'};[{'x0' 'dg' 't0'}' num2cell(round(S,2))]]
    bar(S(:,1:2));
    set(gca,'ylim',[0 6],'xticklabel',{'x0' 'dg' 't0'},'fontsize',16);
    ylabel('mean value');
    title(sprintf('deA > 0, deadenylation (n=%d):\n',sum(k)));
    legend({'ribo' 'polya'},'location','northoutside','box','off');
    saveas(h,'results/param.bar.jpg','jpg');
    saveas(h,'results/param.bar.eps','epsc');

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = [-0.2 0 0.2];
    d = prctile(abs(R(:,4)),20);
    Z(1,:) = [sum(R(:,4)<=-1*d) sum(abs(R(:,4))<d) sum(R(:,4)>=d)];
    Z(2,:) = [sum(R(:,4)<0) sum(abs(R(:,4))==0) sum(R(:,4)>0)];
    a = Z./sum(Z,2);
    bar(a,'stacked');
    L = {};
    S = {'pA' 'NA' 'dA'};
    for i = 1:max(size(x))
        L{i} = sprintf('%s (%.0f%%, %.0f%%)',S{i},100*a(1,i),100*a(2,i));
    end
    set(gca,'xtick',1:2,'xticklabel',{sprintf('T=%.3f',d) sprintf('T=0')},'fontsize',18);
    axis tight;
    legend(L,'location','bestOutside','box','off');
    saveas(h,'results/param.A.jpg','jpg');
    saveas(h,'results/param.A.eps','epsc');

    % fold difference polyA/Total
    if (run_zfish)
        mkdir('results/compare_pa_rz');
        dT = 0.1;
        xT = 0:dT:10;

        mX = zeros(nmtid,max(size(xT)));
        R2r = zeros(nmtid,1);
        for j = 1:nr
            X = [];
            for i = 1:nmtid
                X(i,:) = dg_eval_model_2p(xT,Pr{j}(i,:));
            end
            mX = mX + X;
            R2r = R2r + (isnan(Pr{j}(:,1))==0);
        end
        mX = mX./R2r;
        mX(isnan(mX)) = -10;

        aX = zeros(nmtid,max(size(xT)));
        R2a = zeros(nmtid,1);
        for j = 1:na
            X = [];
            for i = 1:nmtid
                X(i,:) = dg_eval_model_2p(xT,Pa{j}(i,:));
            end
            aX = aX + X;
            R2a = R2a + (isnan(Pa{j}(:,1))==0);
        end
        aX = aX./R2a;
        aX(isnan(aX)) = -10;

        x1 = log2(R(:,1)); % degradation rate
        x1(x1<-2) = -2;
        x2 = aX - mX; % polyA/total ratio
        x2(x2<-3) = -3;
        x2(x2>3) = 3;
        c = corr(x1,x2,'rows','complete');
        num2cell([xT' c'])
        [imx,it] = max(abs(c));
        [imx xT(it)]

        n = size(xT,2);
        for d = 1:100
            w = x2(:,1:(n-d)) - x2(:,1+d:n);
            cd(d) = max(abs(corr(x1,w,'rows','complete')));
        end
        [jmx,jd] = max(cd);
        [jmx dT*jd]

        h = figure;
        scrsz = get(0,'ScreenSize');
        set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
        hold on;
        plot(xT,c,'-k','linewidth',2);
        L = {'A+/total'};
        c2 = corr(x1,x2(:,1:(n-jd))-x2(:,1+jd:n),'rows','complete');
        plot(xT(1:n-jd),c2,'-','linewidth',2);
        L{2} = sprintf('A+/total diff with + %.1f min',60*dT*jd);
        hold off;
        set(gca,'ylim',[-0.35 0.35]);
        xlabel('hpf');
        ylabel('correlation coefficient');
        legend(L,'box','off');
        box('off');
        set(gca,'fontsize',14);
        saveas(h, 'results/compare_pa_rz/dg.pa.1.jpg','jpg');

        clf;
        subplot(2,4,1);
        i = (~isnan(x1)).*(~isnan(x2(:,it))) == 1;
        dscatter(x1(i,1),x2(i,it),'MSIZE',20);
        axis square;
        axis tight;
        xlabel('degradation rate; log2');
        ylabel(sprintf('A+/total ratio at %.2f hr', xT(it)));
        title(sprintf('n=%d, r=%.2f',sum(i),corr(x1(i,1),x2(i,it))));
        set(gca,'fontsize',16);
        subplot(2,4,2);
        [~,ci] = max(c2);
        x3 = x2(:,ci)-x2(:,ci+jd);
        i = (~isnan(x1)).*(~isnan(x3)) == 1;
        dscatter(x1(i,1),x3(i,1),'MSIZE',20);
        axis square;
        axis tight;
        xlabel('degradation rate; log2');
        ylabel(sprintf('A+/total ratio diff (%.1f - %.1f)', xT(ci),xT(ci+jd)));
        title(sprintf('n=%d, r=%.2f',sum(i),corr(x1(i,1),x3(i,1))));
        set(gca,'fontsize',16);
        saveas(h, 'results/compare_pa_rz/dg.pa.2.jpg','jpg');

        % deadenylation rate relative to expression
        h = figure;
        scrsz = get(0,'ScreenSize');
        set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
        c = corr(Rx,R(:,4),'rows','pairwise');
        for i = 1:4
            subplot(2,2,i);
            x = R(:,4);
            y = Rx(:,i);
            y(y<-6) = -6;
            k = (~isnan(x)).*(~isnan(y))==1;
            dscatter(x(k),y(k),'MSIZE',25);
            axis square;
            xlabel('dA rate');
            ylabel(sprintf('polyA/Total ratio (t=%d hr)',i-1));
            title(sprintf('n=%d,r=%.2f',sum(k),c(i)));
            set(gca,'xlim',[-3 3],'ylim',[-6 4],'fontsize',16);
            hold on;
            line([0 0],[-6 4],'LineStyle','-','color','k','linewidth',1);
            hold off;
        end
        saveas(h, 'results/compare_da.jpg','jpg');
    end

    % plot examples
    mkdir('results/model_examples');
    R2r = [];
    for j = 1:nr
        R2r(:,j) = Rr{j};
    end
    R2a = [];
    for j = 1:na
        R2a(:,j) = Ra{j};
    end
    
    minX0 = 4;
    minR2 = min(prctile([R2r R2a],90));%0.8;
    if (nr>0)
        k = (sum(Mr_x0>=minX0,2)==nr);
        k = k.*(sum(R2r>=minR2,2)==nr) == 1;
    else
        k = (sum(Ma_x0>=minX0,2)==na);
    end
    if (na>0)
        k = k.*(sum(R2a>=minR2,2)==na) == 1;
    end
    sum(k)
    if (sum(k)>0)
        j = find(k);
        last = min(sum(k),10);
        plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/fit_',1,0,mzt_fit,1,fpkm_pref);
    end
    
    maxR2 = max(prctile([R2r R2a],10));%0.3;
    if (nr>0)
        k = (sum(Mr_x0>=minX0,2)==nr);
        k = k.*(sum(R2r<maxR2,2)==nr) == 1;
    else
        k = (sum(Ma_x0>=minX0,2)==na);
    end
    if (na>0)
        k = k.*(sum(R2a<maxR2,2)==na) == 1;
    end
    sum(k)
    if (sum(k)>0)
        j = find(k);
        last = min(sum(k),10);
        plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/fitno_',1,0,mzt_fit,0,fpkm_pref);
    end
    
    if (run_frog)
        k = (abs(Pa{1}(:,2))<0.01).*(Ra{1}>=0.95)==1;
        plot_ids(mtid(k),all_data_a,all_data_r,'results/model_examples/test_',1,0,mzt_fit,0,fpkm_pref);
    end

    if (run_zfish)
        k = (sum(Mr_hl>0,2)==nr).*(sum(Mr_hl<0.5,2)==nr).*(sum(Ir,2)==nr)==1;
        sum(k)
        if (sum(k)>0)
            j = find(k);
            last = min(sum(k),10);
            plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/fast_',1,0,mzt_fit,0,fpkm_pref);
        end
        k = (sum(Mr_hl>2,2)==nr).*(sum(Mr_hl<5,2)==nr).*(sum(Ir,2)==nr)==1;
        sum(k)
        if (sum(k)>0)
            j = find(k);
            last = min(sum(k),10);
            plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/slow_',1,0,mzt_fit,0,fpkm_pref);
        end
        k = (sum(Mr_t0>0,2)==nr).*(sum(Mr_t0<2.7,2)==nr).*(sum(Ir,2)==nr)==1;
        sum(k)
        if (sum(k)>0)
            j = find(k);
            last = min(sum(k),10);
            plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/early_',1,0,mzt_fit,0,fpkm_pref);
        end
        k = (sum(Mr_t0>4.5,2)==nr).*(sum(Mr_t0<6,2)==nr).*(sum(Ir,2)==nr)==1;
        sum(k)
        if (sum(k)>0)
            j = find(k);
            last = min(sum(k),10);
            plot_ids(mtid(j(1:last)),all_data_a,all_data_r,'results/model_examples/late_',1,0,mzt_fit,0,fpkm_pref);
        end
    end
end

% ---------------------------------------------------------------------
% expression data: plot data
% ---------------------------------------------------------------------
plot_data = 0;

if (plot_data)
    load('maternal_param.mat','mtid','Pa','Pr');
    
    % expression data
    mxT = 0;
    mT = cell(1,nr);
    mD = cell(1,nr);
    for j = 1:nr
        all_data_r{j}
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_r{j} '.txt']);
        [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        mD{j} = nan(nmtid,size(y,2));
        mD{j}(k1,:) = y;
        mT{j} = x;
        mxT = max(max(x),mxT);
    end
    aT = cell(1,nr);
    aD = cell(1,nr);
    for j = 1:na
        all_data_a{j}
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_a{j} '.txt']);
        [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        aD{j} = nan(nmtid,size(y,2));
        aD{j}(k1,:) = y;
        aT{j} = x;
        mxT = max(max(x),mxT);
    end
    
    I = sum(isnan([mD{:} aD{:}]),2)==0;
    sum(I)
    
    % model estimates
    xT = 0:0.2:mxT;
    mX = zeros(nmtid,max(size(xT)));
    R2r = zeros(nmtid,1);
    for j = 1:nr
        X = [];
        for i = 1:nmtid
            X(i,:) = dg_eval_model_2p(xT, [0 Pr{j}(i,2:5)]);
        end
        mX = mX + X;
        R2r = R2r + (isnan(Pr{j}(:,1))==0);
    end
    mX = mX./R2r;
    mX(isnan(mX)) = -10;
    
    aX = zeros(nmtid,max(size(xT)));
    R2a = zeros(nmtid,1);
    for j = 1:na
        X = [];
        for i = 1:nmtid
            X(i,:) = dg_eval_model_2p(xT, [0 Pa{j}(i,2:5)]);
        end
        aX = aX + X;
        R2a = R2a + (isnan(Pa{j}(:,1))==0);
    end
    aX = aX./R2a;
    aX(isnan(aX)) = -10;
    
    % clustering by model estimates
    [~,i] = cluster_sort([mX aX],2);
    
    % plot: model estimates
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) 0.4*scrsz(3) scrsz(4)]);
    subplot(1,2,1);
    im = imagesc(mX(i,:),[-4 4]);
    im.AlphaData = ~isnan(mX(i,:));
    colormap(gene_colormap(2));
    set(gca,'xtick',1:max(size(xT)),'xticklabel',xT,'ytick',[]);
    xtickangle(90);
    title('total-RNA model');
    subplot(1,2,2);
    im = imagesc(aX(i,:),[-4 4]);
    im.AlphaData = ~isnan(aX(i,:));
    colormap(gene_colormap(2));
    set(gca,'xtick',1:max(size(xT)),'xticklabel',xT,'ytick',[]);
    xtickangle(90);
    title('polyA+ RNA model');
    saveas(h, 'results/expr_all.model.jpg','jpg');
    
    % plot: data, normalized by model X0
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    for j = 1:nr
        subplot(1,max(nr+na+1,8),j);
        x = mD{j} - Pr{j}(:,1);
        x(isnan(x)) = 0;
        im = imagesc(x(i,:),[-4 4]);
        im.AlphaData = ~isnan(x(i,:));
        colormap(gene_colormap(2));
        set(gca,'xtick',1:size(mT{j},2),'xticklabel',mT{j},'ytick',[]);
        xtickangle(90);
        title(regexprep(all_data_r{j},'_',' '));
        set(gca,'color',[0.8 0.8 0.8]);
    end
    for j = 1:na
        subplot(1,max(nr+na+1,8),j+nr);
        x = aD{j} - Pa{j}(:,1);
        x(isnan(x)) = 0;
        im = imagesc(x(i,:),[-4 4]);
        im.AlphaData = ~isnan(x(i,:));
        colormap(gene_colormap(2));
        set(gca,'xtick',1:size(aT{j},2),'xticklabel',aT{j},'ytick',[]);
        xtickangle(90);
        title(regexprep(all_data_a{j},'_',' '));
        set(gca,'color',[0.8 0.8 0.8]);
    end
    subplot(1,max(nr+na+1,8),max(nr+na+1,8));
    im = imagesc(x(i,:),[-4 4]);
    im.AlphaData = ~isnan(x(i,:));
    set(gca,'xtick',[],'ytick',[]);
    colorbar;
    saveas(h, 'results/expr_all.x0.jpg','jpg');

    % plot: data, normalized by time 0
    D = [];
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    for j = 1:nr
        subplot(1,max(nr+na+1,8),j);
        kmin = (mT{j} == min(mT{j}));
        x = mD{j} - nanmean(mD{j}(:,kmin),2);
        im = imagesc(x(i,:),[-4 4]);
        im.AlphaData = ~isnan(x(i,:));
        colormap(gene_colormap(2));
        set(gca,'xtick',1:size(mT{j},2),'xticklabel',mT{j},'ytick',[]);
        xtickangle(90);
        title(regexprep(all_data_r{j},'_',' '));
        set(gca,'color',[0.8 0.8 0.8]);
        D = [D x(i,:)];
    end
    for j = 1:na
        subplot(1,max(nr+na+1,8),j+nr);
        kmin = (aT{j} == min(aT{j}));
        x = aD{j} - nanmean(aD{j}(:,kmin),2);
        im = imagesc(x(i,:),[-4 4]);
        im.AlphaData = ~isnan(x(i,:));
        colormap(gene_colormap(2));
        set(gca,'xtick',1:size(aT{j},2),'xticklabel',aT{j},'ytick',[]);
        xtickangle(90);
        title(regexprep(all_data_a{j},'_',' '));
        set(gca,'color',[0.8 0.8 0.8]);
        D = [D x(i,:)];
    end
    subplot(1,max(nr+na+1,8),max(nr+na+1,8));
    im = imagesc(x(i,:),[-4 4]);
    im.AlphaData = ~isnan(x(i,:));
    set(gca,'xtick',[],'ytick',1:100:nmtid);
    colorbar;
    saveas(h, 'results/expr_all.t0.jpg','jpg');

    [~,w] = ismember(mtid(i),tid); % [tid(w(w>0)) mtid(i)]
    write_text_file('results/expr_all.t0.txt',[mtid(i) gid(w(w>0)) num2cell(D)]);

    close all;
end

% ---------------------------------------------------------------------
% fold change
% ---------------------------------------------------------------------
plot_fold = 0;

if (plot_fold)
    load('maternal_param.avg.mat','mtid');
    
    mT = cell(1,nr);
    mD = cell(1,nr);
    mA = cell(1,nr);
    aT = cell(1,nr);
    aD = cell(1,nr);
    aA = cell(1,nr);
    for j = 1:nr
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_r{j} '.txt']);
        [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        mD{j} = nan(nmtid,size(y,2));
        mD{j}(k1,:) = y;
        mT{j} = x;
        [~,~,u] = unique(mT{j});
        for k = 1:size(mD{j},1)
            mA{j}(k,:) = accumarray(u,mD{j}(k,:),[],@nanmean);
        end
    end
    for j = 1:na
        [D,t,iM,iP,gid,tid,gc] = load_data([fpkm_pref all_data_a{j} '.txt']);
        [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
        x = t(iM);
        y = D(k2(k2>0),iM);
        aD{j} = nan(nmtid,size(y,2));
        aD{j}(k1,:) = y;
        aT{j} = x;
        [~,~,u] = unique(aT{j});
        for k = 1:size(aD{j},1)
            aA{j}(k,:) = accumarray(u,aD{j}(k,:),[],@nanmean);
        end
    end
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = -10:0.5:10;
    hold on;
    mF = [];
    for j = 1:nr
        t = unique(mT{j});
        kmax = (t == max(t));
        kmin = (t == min(t));
        dt = max(t) - min(t);
        dx = mA{j}(:,kmax) - mA{j}(:,kmin);
        mF(:,j) = dx;
    end
    y = hist(mF,x);
    plot(x,y./sum(y),'-','linewidth',2);
    aF = [];
    for j = 1:na
        t = unique(aT{j});
        kmax = (t == max(t));
        kmin = (t == min(t));
        dt = max(t) - min(t);
        dx = aA{j}(:,kmax) - aA{j}(:,kmin);
        aF(:,j) = dx;
    end
    y = hist(aF,x);
    plot(x,y./sum(y),':','linewidth',2);
    line([0 0],[0 0.15],'LineStyle','-','color','k','linewidth',1);
    hold off;
    xlabel('fold-change: t(max) - t(min)');
    ylabel('fraction of genes');
    legend(regexprep([all_data_r all_data_a],'_',' '),'box','off','location','bestOutside');
    set(gca,'fontsize',18);
    saveas(h, 'results/expr_fold_change.1.jpg','jpg');
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = 0:0.5:10;
    hold on;
    mF = [];
    for j = 1:nr
        mF(:,j) = max(mA{j},[],2) - min(mA{j},[],2);
    end
    y = hist(mF,x);
    plot(x,y./sum(y),'-','linewidth',2);
    aF = [];
    for j = 1:na
        aF(:,j) = max(aA{j},[],2) - min(aA{j},[],2);
    end
    y = hist(aF,x);
    plot(x,y./sum(y),':','linewidth',2);
    line([0 0],[0 0.15],'LineStyle','-','color','k','linewidth',1);
    hold off;
    xlabel('fold-change: max - min');
    ylabel('fraction of genes');
    legend(regexprep([all_data_r all_data_a],'_',' '),'box','off','location','bestOutside');
    set(gca,'fontsize',18);
    saveas(h, 'results/expr_fold_change.2.jpg','jpg');
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = 0:10;
    for j = 1:nr
        subplot(2,max(maxN,3),j);
        u = unique(mT{j});
        [~,t] = max(mA{j},[],2);
        y = hist(u(t)',x);
        bar(x,y./sum(y));
        title(regexprep(all_data_r(j),'_',' '));
        xlabel('time of maximal expression');
        ylabel('fraction of genes');
        set(gca,'ylim',[0 1],'xlim',[-1 11],'xtick',0:10,'fontsize',12);
    end
    for j = 1:na
        subplot(2,max(maxN,3),nr+j);
        u = unique(aT{j});
        [~,t] = max(aA{j},[],2);
        y = hist(u(t)',x);
        bar(x,y./sum(y));
        title(regexprep(all_data_a(j),'_',' '));
        xlabel('time of maximal expression');
        ylabel('fraction of genes');
        set(gca,'ylim',[0 1],'xlim',[-1 11],'xtick',0:10,'fontsize',12);
    end
    saveas(h, 'results/expr_fold_change.3.jpg','jpg');

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = -4:4;
    hold on;
    N = zeros(2,3);
    y = [];
    k = 1;
    for j = 1:nr
        t = unique(mT{j});
        kmin = (t == min(t));
        dt = max(t) - min(t);
        dx = mD{j}(:,kmin==0) - mA{j}(:,kmin);
        y(:,k) = hist(dx(:),x);
        k = k+1;
        N(1,1) = N(1,1)+sum(abs(dx(:))<log2(1.2));
        N(1,2) = N(1,2)+sum(dx(:)<-1*log2(1.2));
        N(1,3) = N(1,3)+sum(dx(:)>log2(1.2));
    end
    for j = 1:na
        t = unique(aT{j});
        kmin = (t == min(t));
        dt = max(t) - min(t);
        dx = aD{j}(:,kmin==0) - aA{j}(:,kmin);
        y(:,k) = hist(dx(:),x);
        k = k+1;
        N(2,1) = N(2,1)+sum(abs(dx(:))<log2(1.2));
        N(2,2) = N(2,2)+sum(dx(:)<-1*log2(1.2));
        N(2,3) = N(2,3)+sum(dx(:)>log2(1.2));
    end
    b = bar(x,y./sum(y));
    for j = 1:nr
        b(j).FaceColor = [.8 .8 .8];
    end
    for j = 1:na
        b(j+nr).FaceColor = [.2 .2 .9];
    end
    line([0 0],[0 1],'LineStyle','-','color','k','linewidth',1);
    hold off;
    xlabel('fold-change relative to time 0');
    ylabel('fraction of data');
    legend(regexprep([all_data_r all_data_a],'_',' '),'box','off','location','bestOutside');
    set(gca,'fontsize',18);
    axis tight;
    saveas(h, 'results/expr_fold_change.4.jpg','jpg');

    clf;
    c = N./sum(N,2);
    bar((c)');
    for i = 1:3
        text(i-0.2,c(1,i),sprintf('%.0f %%',100*c(1,i)),'fontsize',18);
        text(i+0.2,c(2,i),sprintf('%.0f %%',100*c(2,i)),'fontsize',18);
    end
    set(gca,'xtick',1:3,'xticklabel',{'less than 1.2 fold' 'down > 1.2 fold' 'up > 1.2 fold'});
    set(gca,'ylim',[0 1],'fontsize',18);
    ylabel('fraction of data');
    legend({'total-RNA' 'polyA+RNA'},'box','off','location','bestOutside');
    saveas(h, 'results/expr_fold_change.4all.jpg','jpg');
    saveas(h, 'results/expr_fold_change.4all.eps','epsc');
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    x = -3:0.5:7;
    subplot(1,2,1);
    hold on;
    T = [];
    mF0 = [];
    for j = 1:nr
        t = unique(mT{j});
        kmax = (t == max(t));
        mF0(:,j) = mA{j}(:,kmax);
        T = [T max(t)];
    end
    y = hist(mF0,x);
    plot(x,y./sum(y),'-','linewidth',2);
    aF0 = [];
    for j = 1:na
        t = unique(aT{j});
        kmax = (t == max(t));
        aF0(:,j) = aA{j}(:,kmax);
        T = [T max(t)];
    end
    y = hist(aF0,x);
    plot(x,y./sum(y),':','linewidth',2);
    hold off;
    xlabel('late time expression');
    ylabel('fraction of genes');
    L = regexprep(strcat(regexprep([all_data_r all_data_a],'_',' ')',' (',num2str(T'),' hr)'),'(  *','(');
    legend(L,'box','off','location','northOutside');
    set(gca,'fontsize',18);
    axis tight;
    subplot(1,2,2);
    hold on;
    T = [];
    mF1 = [];
    for j = 1:nr
        t = unique(mT{j});
        kmin = (t == min(t));
        mF1(:,j) = mA{j}(:,kmin);
        T = [T min(t)];
    end
    y = hist(mF1,x);
    plot(x,y./sum(y),'-','linewidth',2);
    aF1 = [];
    for j = 1:na
        t = unique(aT{j});
        kmin = (t == min(t));
        aF1(:,j) = aA{j}(:,kmin);
        T = [T min(t)];
    end
    y = hist(aF1,x);
    plot(x,y./sum(y),':','linewidth',2);
    hold off;
    xlabel('initial time expression');
    ylabel('fraction of genes');
    L = regexprep(strcat(regexprep([all_data_r all_data_a],'_',' ')',' (',num2str(T'),' hr)'),'(  *','(');
    legend(L,'box','off','location','northOutside');
    set(gca,'fontsize',18);
    axis tight;
    saveas(h, 'results/expr_fold_change.5.jpg','jpg');
    
    clf;
    subplot(1,2,1);
    hold on;
    y = hist(mF0(:),x);
    plot(x,y./sum(y),'-','linewidth',2);
    y = hist(aF0(:),x);
    plot(x,y./sum(y),'-','linewidth',2);
    hold off;    
    xlabel('late time expression');
    ylabel('fraction of genes');
    L = {sprintf('total RNA (%.2f)',2.^nanmean(mF0(:))) sprintf('polyA+ RNA (%.2f)',nanmean(aF0(:)))};
    legend(L,'box','off','location','northOutside');
    title(sprintf('fold = %.2f',2.^[nanmean(mF0(:))-nanmean(aF0(:))]));
    set(gca,'fontsize',18);
    axis tight;
    subplot(1,2,2);
    hold on;
    y = hist(mF1(:),x);
    plot(x,y./sum(y),'-','linewidth',2);
    y = hist(aF1(:),x);
    plot(x,y./sum(y),'-','linewidth',2);
    hold off;
    xlabel('initial time expression');
    ylabel('fraction of genes');
    L = {sprintf('total RNA (%.2f)',2.^nanmean(mF1(:))) sprintf('polyA+ RNA (%.2f)',nanmean(aF1(:)))};
    legend(L,'box','off','location','northOutside');
    title(sprintf('fold = %.2f',2.^[nanmean(mF1(:))-nanmean(aF1(:))]));
    set(gca,'fontsize',18);
    axis tight;
    saveas(h, 'results/expr_fold_change.5all.jpg','jpg');
    saveas(h, 'results/expr_fold_change.5all.eps','epsc');

    close all;
end

% ---------------------------------------------------------------------
% plot known kmers
% ---------------------------------------------------------------------
plot_kmer_known = 0;
plot_kmer_known_individual_param = 0;

if (plot_kmer_known)
    Kids = {'GCACTT' 'TATTTAT' 'TATTAT' 'TGTA.ATA' ...
        'TTTTT' 'TTTTTT' 'TTTTTTT' ...
        'TTTTAT' 'TTTTTA' 'TTTTAGT' ...
        'AAAA' 'AAAAA' 'AAAAAA' 'AAAAAAA' ...
        'AATAAA' 'AAAATA' 'AAATAA' 'ATAAAA'...
        'AAGAAA' 'AAAAGA' 'AAAGAA' 'AGAAAA' ...
        'CCCC' 'CCCCC' 'CCCCCC' ...
        'GGGG' 'GGGGG' 'GGGGGG' ...
        'CTCC|CCTC' 'CACA|ACAC' ...
        'CCCCG' 'CCGGG' 'CCCCAG' 'CTGC|CCTG|CTGG|GCTG'}; % xenopus
    load('maternal_param.avg.mat','mtid','R','Rx','Rx_id', ...
        'Mr_dg','Ma_dg','Ma_da','Mr_t0','Ma_t0','Mr_x0','Ma_x0','Xa','Xr');

    % load sequences
    f = fopen(seq_file);
    X = textscan(f,'%s %s');
    fclose(f);
    S = X{2};
    Sid = regexprep(X{1},'.*;','');
    
    [~,j1,j2] = intersect(mtid,Sid);
    Sid = Sid(j2(j2>0));
    S = S(j2(j2>0));

    % for i = 1:max(size(Kids))
    %     k = cellfun(@isempty,regexp(S,Kids{i})) == 0;
    %     write_text_file(['results/kmer_known/' Kids{i} '_genes.txt'], Sid(k));
    % end

    % all
    mkdir('results/kmer_known');
    Di = {'x0r' 'x0a' 't0r' 't0a' 'dgr' 'dga' 'da'};
    D = [nanmean(Mr_x0,2) nanmean(Ma_x0,2) nanmean(Mr_t0,2) nanmean(Ma_t0,2) nanmean(Mr_dg,2) nanmean(Ma_dg,2) nanmean(Ma_da,2)];
    D = scale_param(D);
    k = sum(~isnan(D)) > 10;
    plot_known_kmers(Sid,S,Kids,D(j1,k),Di(k),'results/kmer_known/all','scaled param (a.u.)',[0 10],0.5);
    plot_known_kmers(Sid,S,Kids,D(j1,k),Di(k),'results/kmer_known/all.box','scaled param (a.u.)',[0 10],0.5,4,1);

    % Ratio
    if ((nr>0)*(na>0)==1)
        plot_known_kmers(Sid,S,Kids,Rx(j1,:),Rx_id,'results/kmer_known/r','A/T ratio',[-2 1],0.5);
        plot_known_kmers(Sid,S,Kids,Rx(j1,:),Rx_id,'results/kmer_known/r','A/T ratio',[-2 1],0.5,4,1);
    end

    if ((nr+na > 1) && (plot_kmer_known_individual_param))
        % X0l
        D = [Mr_x0 Ma_x0];
        Di = [all_data_r all_data_a];
        plot_known_kmers(Sid,S,Kids,D(j1,1:(nr+na)),Di,'results/kmer_known/x0','logX0',[-5 12],0.5);
        
        % dA rate
        if (na > 0)
            D = Ma_da;
            Di = all_data_a;
            plot_known_kmers(Sid,S,Kids,D(j1,1:na),Di,'results/kmer_known/da','deA rate (1/hr)',[-3 3],0.2);
        end

        % degradation rate
        D = [Mr_dg Ma_dg];
        Di = [all_data_r all_data_a];
        plot_known_kmers(Sid,S,Kids,D(j1,1:(nr+na)),Di,'results/kmer_known/dg','dg rate (1/hr)',[0 3],0.1);
        
        % onset time
        D = [Mr_t0 Ma_t0];
        Di = [all_data_r all_data_a];
        plot_known_kmers(Sid,S,Kids,D(j1,1:(nr+na)),Di,'results/kmer_known/t0','onset (hr)',[0 7],0.5);
    end
end

% ---------------------------------------------------------------------
% expression data: kmer enrichments
% ---------------------------------------------------------------------
plot_kmer = 1;
plot_kmer_individual_param = 0;

if (plot_kmer)
    load('maternal_param.avg.mat','mtid','R','Rx','Rx_id', ...
        'Mr_dg','Ma_dg','Ma_da','Mr_t0','Ma_t0','Mr_x0','Ma_x0');
    
    % load sequences
    f = fopen(seq_file);
    X = textscan(f,'%s %s');
    fclose(f);
    S = X{2};
    Sid = regexprep(X{1},'.*;','');
    
    [~,j1,j2] = intersect(mtid,Sid); % [mtid(j1) Sid(j2(j2>0))]
    Sid = Sid(j2(j2>0));
    S = S(j2(j2>0));
    
    % all
    Pid = {'x0r' 'x0a' 't0r' 't0a' 'dgr' 'dga' 'da'};
    Pall = [nanmean(Mr_x0,2) nanmean(Ma_x0,2) nanmean(Mr_t0,2) nanmean(Ma_t0,2) nanmean(Mr_dg,2) nanmean(Ma_dg,2) nanmean(Ma_da,2)];
    k = sum(~isnan(Pall)) > 10;
    kmers_all(mtid(j1),S,Pall(j1,k),Pid(k),'results/kmer_all',kmer_range,kmer_dir,kmer_esize,kmer_alpha);

    % Ratio
    if ((nr>0)*(na>0)==1)
        kmers_all(mtid(j1),S,Rx(j1,:),Rx_id,'results/kmer_r',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
    end

    if ((nr>1 || na>1) && (plot_kmer_individual_param))
        % X0
        if (na>0)
            x = [Mr_x0 Ma_x0];
            kmers_all(mtid(j1),S,x(j1,:),[all_data_r all_data_a],'results/kmer_x0',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        else
            kmers_all(mtid(j1),S,Mr_x0(j1,:),all_data_r,'results/kmer_x0',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        end
        
        % pA rate
        if (na>0)
            pa = Ma_da;
            kmers_all(mtid(j1),S,pa(j1,:),all_data_a,'results/kmer_da',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        end
        
        % degradation rate
        if (na>0)
            d = [Mr_dg Ma_dg];
            kmers_all(mtid(j1),S,d(j1,:),[all_data_r all_data_a],'results/kmer_dg',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        else
            kmers_all(mtid(j1),S,Mr_dg(j1,:),all_data_r,'results/kmer_dg',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        end
        
        % onset time
        % less accurate so require a larger effect size
        if (na>0)
            t = [Mr_t0 Ma_t0];
            kmers_all(mtid(j1),S,t(j1,:),[all_data_r all_data_a],'results/kmer_t0',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        else
            kmers_all(mtid(j1),S,Mr_t0(j1,:),all_data_r,'results/kmer_t0',kmer_range,kmer_dir,kmer_esize,kmer_alpha);
        end   
    end
end


function run_fit(mtid,all_data,mzt,pref)

mN = size(mtid,1);
L = {'ensid' 'x0' 'dg' 't0' 't1' 'rsq' 'mse'};
for j = 1:max(size(all_data))
    [D,t,iM,~,~,tid] = load_data([pref all_data{j} '.txt']);    
    [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = D(k2(k2>0),iM);
    
    fprintf('fit model 1p: %s\n', all_data{j});
    if (exist(['model_param.1p.' all_data{j} '.mat'],'file'))
        load(['model_param.1p.' all_data{j} '.mat'],'mtid','mP','mR','mE');
    else
        P = nan(size(y,1),4);
        E = nan(size(y,1),1);
        R = nan(size(y,1),1);
        if (size(x,2)>=4)
            parfor i = 1:size(y,1)
                [P(i,:),E(i,:),~,R(i,:)] = fit_model(x,y(i,:),mzt);
            end
        end
        mP = nan(mN,size(P,2));
        mP(k1,:) = P;
        mR = nan(mN,size(R,2));
        mR(k1,:) = R;
        mE = nan(mN,size(E,2));
        mE(k1,:) = E;
        save(['model_param.1p.' all_data{j} '.mat'],'mtid','mP','mR','mE');
        write_text_file(['model_param.1p.' all_data{j} '.txt'],[L;[mtid num2cell([mP mR mE])]]);
    end
    k = sum(~isnan(mP(:,1)));
    fprintf('fitted model 1p: fit %d out of %d (%.1f%%)\n',k,mN,100*k/mN);
end    

function [mP,mR,mE] = run_fit_2p(mtid,all_data,mzt,pref)

mN = size(mtid,1);
L = {'ensid' 'x0' 'da' 'dg' 't0' 't1' 'rsq' 'mse'};
for j = 1:max(size(all_data))
    [D,t,iM,~,~,tid] = load_data([pref all_data{j} '.txt']);
    [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = D(k2(k2>0),iM);
    
    fprintf('fit model 2p: %s\n', all_data{j});
    if (exist(['model_param.2p.' all_data{j} '.mat'],'file'))
        load(['model_param.2p.' all_data{j} '.mat'],'mtid','mP','mR','mE');
    else
        P = nan(size(y,1),5);
        E = nan(size(y,1),1);
        R = nan(size(y,1),1);
        if (size(x,2)>=5)
            parfor i = 1:size(y,1)
                [P(i,:),E(i,:),~,R(i,:)] = fit_model_2p(x,y(i,:),mzt);
            end
        end
        mP = nan(mN,size(P,2));
        mP(k1,:) = P;
        mR = nan(mN,size(R,2));
        mR(k1,:) = R;
        mE = nan(mN,size(E,2));
        mE(k1,:) = E;
        save(['model_param.2p.' all_data{j} '.mat'],'mtid','mP','mR','mE');
        write_text_file(['model_param.2p.' all_data{j} '.txt'],[L;[mtid num2cell([mP mR mE])]]);
    end
    k = sum(~isnan(mP(:,1)));
    fprintf('fitted model 2p: fit %d out of %d (%.1f%%)\n',k,mN,100*k/mN);
end

function nD = scale_param(D)
% D = ['x0r' 'x0a' 't0r' 't0a' 'dgr' 'dga' 'da']

xmax = max(max(D(:,[1 2])));
xmin = min(min(D(:,[1 2])));
tmax = max(max(D(:,[3 4])));
tmin = 0;
dmax = max(max(D(:,[5 6])));
dmin = min(min(D(:,[5 6])));
rmax = max(abs(D(:,7)));
rmin = -1*rmax;
Dmin = [xmin xmin tmin tmin dmin dmin rmin];
Dmax = [xmax xmax tmax tmax dmax dmax rmax];
nD = 10*(D - Dmin)./(Dmax - Dmin);
