function kmers_all(Dids,Dseq,Dval,Dnames,output_dir,Krange,Kdir,esize,palpha,erange,prange,sort_type,suffix)

if (nargin < 6)
    Krange = 3:7;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    esize = 0.3;
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
if (nargin < 13)
    suffix = '';
end

n = size(Dval,2);
minP = 1e-200;
min_overlap = min(Krange);
mxP = prange(2);
mxE = max(abs(erange));

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

% KS test to identify significant kmers
% W1/W2 = [kmer] [pvalue neg. E] [pvalue pos. E] [gene +kmer count] [effect size]
if (exist([output_dir '/ks.mat'],'file'))
    load([output_dir '/ks.mat'],'W1','W2');
    fprintf('load KS test data\n');
else
    W1 = cell(1,n);
    W2 = cell(1,n);
    for i = 1:n
        Dnames{i}
        I = (~isnan(Dval(:,i)));
        [W1{i},W2{i}] = volcano_ks(Dids(I),Dval(I,i),[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
    end
    save([output_dir '/ks.mat'],'W1','W2');
end

% filter: 1. kmers that are contained within another kmer with a more significant p-value
%         2. kmers that contain another kmer with a more significant p-value
% P1/2 = [kmer] [pvalue] [effect size] [no. of kmers merged]
if (exist([output_dir '/kmer.mat'],'file'))
    load([output_dir '/kmer.mat'],'P1','P2');
    fprintf('load KMER data\n');
else
    P1 = cell(1,n);
    P2 = cell(1,n);
    for i = 1:n
        if (~isempty(W1{i}))
            P1{i} = filter_kmers_by_seq(W1{i}(:,1),cell2mat(W1{i}(:,[2 5])),min_overlap);
        end
        if (~isempty(W2{i}))
            P2{i} = filter_kmers_by_seq(W2{i}(:,1),cell2mat(W2{i}(:,[3 5])),min_overlap);
        end
    end
    save([output_dir '/kmer.mat'],'P1','P2');
end

% combine kmers from all different datasets
% ** Use only the identity of filtered k-mers, not p-values 
%    (take pvalues from the original set)
if (exist([output_dir '/kmer_all.mat'],'file'))
    load([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all');
    fprintf('load KMER combined\n');
else
    S_all = [];
    for i = 1:n
        if (~isempty(P1{i}))
            S_all = [S_all;P1{i}(:,1)];
        end
        if (~isempty(P2{i}))
            S_all = [S_all;P2{i}(:,1)];
        end
    end
    S_all = unique(S_all);
    Q_all = nan(size(S_all,1),n);
    E_all = nan(size(S_all,1),n);
    for i = 1:n
        if (~isempty(P1{i}))
            [w1,w2] = ismember(S_all,W1{i}(:,1)); % [all_seq(w1) W1{i}(w2(w2>0),1)]
            if (sum(w1)>0)
                Q_all(w1,i) = cell2mat(W1{i}(w2(w2>0),2));
                E_all(w1,i) = cell2mat(W1{i}(w2(w2>0),5));
            end
        end
        if (~isempty(P2{i}))
            [w1,w2] = ismember(S_all,W2{i}(:,1)); % [all_seq(w1) W1{i}(w2(w2>0),1)]
            if (sum(w1)>0)
                Q_all(w1,i) = cell2mat(W2{i}(w2(w2>0),3));
                E_all(w1,i) = cell2mat(W2{i}(w2(w2>0),5));
            end
        end
    end

    if (~isempty(S_all))
        P_all = -1*log10(Q_all);
        P_all(P_all > -1*log10(minP)) = -1*log10(minP);

        % volcano plots
        h = figure;
        scrsz = get(0,'ScreenSize');
        set(h, 'OuterPosition',[1 scrsz(4) 0.5*scrsz(3) scrsz(4)]);
        for j = 1:n
            clf;
            k0 = isnan(P_all(:,j)) + isnan(E_all(:,j)) == 0;
            k1 = (P_all(:,j) >= -1*log10(palpha)).*(abs(E_all(:,j)) >= esize) == 1;
            mxPj = max(max(P_all(:,j)),-1*log10(palpha)) + 1;
            mxEj = max(max(abs(E_all(:,j))),esize) + 0.1;
            hold on;
            for i = find(k1)'
                seq = regexprep(S_all(i),'T','U');
                text(E_all(i,j),P_all(i,j),seq,'fontsize',14);
            end
            if (sum(k0) > 10)
                n1 = zeros(size(E_all(k0,j)));
                n2 = zeros(size(P_all(k0,j)));
                if (max(size(unique(E_all(k0,j))))<2)
                    n1 = 0.01*randn(size(E_all(k0,j)));
                end
                if (max(size(unique(P_all(k0,j))))<2)
                    n2 = 0.01*randn(size(P_all(k0,j)));
                end
                dscatter(E_all(k0,j)+n1,P_all(k0,j)+n2,'MSIZE',50);
            else
                plot(E_all(k0,j),P_all(k0,j),'.','markerSize',20);
            end
            line([0 0],[0 mxPj],'LineStyle','-','color','k','linewidth',2);
            line([-1*mxEj mxEj],-1*log10([palpha palpha]),'LineStyle','-','color','r','linewidth',1);
            line([esize esize],[0 mxPj],'LineStyle','-','color','r','linewidth',1);
            line(-1*[esize esize],[0 mxPj],'LineStyle','-','color','r','linewidth',1);
            hold off;
            axis tight;
            axis square;
            box on;
            xlabel('effect size (SMD)');
            ylabel('p-value');
            title(sprintf('%s (n=%d,n=%d)',Dnames{j},sum(k0),sum(k1)));
            set(gca,'xlim',[-1*mxEj mxEj],'ylim',[0 mxPj]);
            saveas(h, [output_dir '/matrix.' Dnames{j} '.volcano.jpg'],'jpg');
            saveas(h, [output_dir '/matrix.' Dnames{j} '.volcano.eps'],'epsc');

            if (sum(k0)>0)
                W = sortrows([S_all(k0) num2cell([10.^(-1*P_all(k0,j)) E_all(k0,j)])],[2 3]);
                write_text_file([output_dir '/matrix.' Dnames{j} '.volcano.txt'],W);
            end
            if (sum(k1)>0)
                W = sortrows([S_all(k1) num2cell([10.^(-1*P_all(k1,j)) E_all(k1,j)])],[2 3]);
                write_text_file([output_dir '/matrix.' Dnames{j} '.volcano.e.txt'],W);
            end
        end
        close all;
    else
        P_all = [];
    end

    save([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all');
end
fprintf('total: %d kmers (all datasets)\n', size(S_all,1));

% plot results
if (~isempty(S_all))
    P = P_all;
    P(E_all<0) = -1*P_all(E_all<0);
    P(isnan(P_all)) = 0;
    E = E_all;
    E(isnan(E_all)) = 0;

    % all kmers that pass both pvalue and esize thresholds
    k = zeros(size(P,1),1);
    for j = 1:n
        pj = (abs(P(:,j))>=-1*log10(palpha)).*(abs(E(:,j))>=esize) == 1;
        k = k + pj > 0;
    end
    Sk = S_all(k,:);
    Ek = E(k,:);
    Pk = P(k,:);
    if (sum(k) <= 500)
        [h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),esize,mxP,mxE);
        set(gca,'xtick',[]);
        saveas(h, [output_dir '/matrix.all.jpg'],'jpg');
        write_text_file([output_dir '/matrix.P.all.txt'],[['kmer' Dnames];[Sk(i) num2cell(Pk(i,:))]]);
        write_text_file([output_dir '/matrix.E.all.txt'],[['kmer' Dnames];[Sk(i) num2cell(Ek(i,:))]]);
    else
        write_text_file([output_dir '/matrix.P.all.txt'],[['kmer' Dnames];[Sk num2cell(Pk)]]);
        write_text_file([output_dir '/matrix.E.all.txt'],[['kmer' Dnames];[Sk num2cell(Ek)]]);
    end
    
    % top N kmers per column
    plot_Nkmers(20,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir);
    plot_Nkmers(10,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir);
    
    % top 50 kmers per column (separate)
    for j = 1:n
        plot_Nkmers_col(j,50,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir);
    end
    close all;
end


function plot_Nkmers(Nkmer,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir)

n = size(P,2);
k = zeros(size(P,1),1)==1;
for j = 1:n
    pj = abs(P(:,j));
    ej = abs(E(:,j));
    i = (pj>=-1*log10(palpha)).*(ej>=esize) == 1;
    pi = sort(pj(i));
    if (size(pi,1)>=Nkmer)
        k = k + (pj>=min(pi(end-Nkmer+1:end))).*(ej>=esize) > 0;
    elseif (~isempty(pi))
        k = (k + i) > 0;
    end
end
Sk = S_all(k,:);
Ek = E(k,:);
Pk = P(k,:);

[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE);
set(gca,'xtick',[]);
saveas(h, [output_dir '/matrix.N' num2str(Nkmer) '.jpg'],'jpg');
write_text_file([output_dir '/matrix.P.N' num2str(Nkmer) '.txt'],[Sk(i) num2cell(Pk(i,:))]);
write_text_file([output_dir '/matrix.E.N' num2str(Nkmer) '.txt'],[Sk(i) num2cell(Ek(i,:))]);


function plot_Nkmers_col(j,Nkmer,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir)

pj = abs(P(:,j));
ej = abs(E(:,j));
i = (pj>=-1*log10(palpha)).*(ej>=esize) == 1;
pi = sort(pj(i));
if (size(pi,1)>=Nkmer)
    k = (pj>=min(pi(end-Nkmer+1:end))).*(ej>=esize) > 0;
elseif (~isempty(pi))
    k = i > 0;
else
    k = zeros(size(P,1),1)==1;
end
Sk = S_all(k,:);
Ek = E(k,:);
Pk = P(k,:);

[~,I] = sort(Pk(:,j));
[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE,I);
set(gca,'xtick',[]);
saveas(h, [output_dir '/matrix.' Dnames{j} '.jpg'],'jpg');
write_text_file([output_dir '/matrix.P.' Dnames{j} '.txt'],[Sk(i) num2cell(Pk(i,:))]);
write_text_file([output_dir '/matrix.E.' Dnames{j} '.txt'],[Sk(i) num2cell(Ek(i,:))]);
