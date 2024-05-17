function kmers_all_merge(dir_list,output_dir,prefix,esize,palpha,erange,prange,order,sort_type)

if (nargin < 4)
    esize = 0.3;
end
if (nargin < 5)
    palpha = 0.01;
end
if (nargin < 6)
    erange = [-1 1];
end
if (nargin < 7)
    prange = [0 30];
end
if (nargin < 8)
    order = [];
end
if (nargin < 9)
    sort_type = 1; % 1 = by pvalues; 2 = by kmers
end

minP = 1e-200;
mxP = prange(2);
mxE = max(abs(erange));

mkdir(output_dir);


% combine kmers from all different datasets
% ** Use only the identity of filtered k-mers, not p-values 
%    (take pvalues from the original set)
if (exist([output_dir '/kmer_all.mat'],'file'))
    load([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all','n');
    fprintf('load KMER combined\n');
else
    % collect all kmers
    S_all = [];
    n = [];
    for d = 1:max(size(dir_list))
        load([dir_list{d} '/kmer.mat'],'P1','P2');
        n(d) = size(P1,2);
        for i = 1:n(d)
            if (~isempty(P1{i}))
                S_all = [S_all;P1{i}(:,1)];
            end
            if (~isempty(P2{i}))
                S_all = [S_all;P2{i}(:,1)];
            end
        end
    end
    S_all = unique(S_all);

    % collect all pvalues
    Q_all = nan(size(S_all,1),sum(n));
    E_all = nan(size(S_all,1),sum(n));
    k = 1;
    for d = 1:max(size(dir_list))
        load([dir_list{d} '/kmer.mat'],'P1','P2');
        load([dir_list{d} '/ks.mat'],'W1','W2');
        for i = 1:n(d)
            if (~isempty(P1{i}))
                [w1,w2] = ismember(S_all,W1{i}(:,1)); % [all_seq(w1) W1{i}(w2(w2>0),1)]
                if (sum(w1)>0)
                    Q_all(w1,k) = cell2mat(W1{i}(w2(w2>0),2));
                    E_all(w1,k) = cell2mat(W1{i}(w2(w2>0),5));
                end
            end
            if (~isempty(P2{i}))
                [w1,w2] = ismember(S_all,W2{i}(:,1)); % [all_seq(w1) W1{i}(w2(w2>0),1)]
                if (sum(w1)>0)
                    Q_all(w1,k) = cell2mat(W2{i}(w2(w2>0),3));
                    E_all(w1,k) = cell2mat(W2{i}(w2(w2>0),5));
                end
            end
            k = k+1;
        end
    end

    if (~isempty(S_all))
        P_all = -1*log10(Q_all);
        P_all(P_all > -1*log10(minP)) = -1*log10(minP);
    else
        P_all = [];
    end
    save([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all','n');
end

fprintf('input: %d kmers, %d properties\n',size(P_all));
Dnames = regexprep(cellstr(num2str((1:sum(n))'))',' ','');
if (~isempty(order))
    P_all = P_all(:,order);
    E_all = E_all(:,order);
    Dnames = Dnames(:,order);

    k = sum(isnan(P_all),2)<size(P_all,2);
    P_all = P_all(k,:);
    E_all = E_all(k,:);
    S_all = S_all(k,:);
end
fprintf('filtered: %d kmers, %d properties\n',size(P_all));

% plot results
if (~isempty(S_all))
    P = P_all;
    P(E_all<0) = -1*P_all(E_all<0);
    P(isnan(P_all)) = 0;
    E = E_all;
    E(isnan(E_all)) = 0;

    % all kmers that pass both pvalue and esize thresholds
    k = zeros(size(P,1),1);
    for j = 1:size(P,2)
        pj = (abs(P(:,j))>=-1*log10(palpha)).*(abs(E(:,j))>=esize) == 1;
        k = k + pj > 0;
    end
    Sk = S_all(k,:);
    Ek = E(k,:);
    Pk = P(k,:);
    if (sum(k) <= 500)
        [h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),esize,mxP,mxE);
        set(gca,'xtick',[]);
        saveas(h, [output_dir '/matrix.' prefix '.all.jpg'],'jpg');
        write_text_file([output_dir '/matrix.' prefix '.P.all.txt'],[['kmer' Dnames];[Sk(i) num2cell(Pk(i,:))]]);
        write_text_file([output_dir '/matrix.' prefix '.E.all.txt'],[['kmer' Dnames];[Sk(i) num2cell(Ek(i,:))]]);
    else
        write_text_file([output_dir '/matrix.' prefix '.P.all.txt'],[['kmer' Dnames];[Sk num2cell(Pk)]]);
        write_text_file([output_dir '/matrix.' prefix '.E.all.txt'],[['kmer' Dnames];[Sk num2cell(Ek)]]);
    end
    
    % top N kmers per column
    plot_Nkmers(20,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir,prefix);
    plot_Nkmers(10,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir,prefix);
    plot_Nkmers(5,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir,prefix);
    
    close all;
end


function plot_Nkmers(Nkmer,P,E,S_all,Dnames,sort_type,palpha,esize,mxP,mxE,output_dir,prefix)

n = size(P,2);
k = zeros(size(P,1),1);
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
saveas(h, [output_dir '/matrix.' prefix '.N' num2str(Nkmer) '.jpg'],'jpg');
write_text_file([output_dir '/matrix.' prefix '.P.N' num2str(Nkmer) '.txt'],[Sk(i) num2cell(Pk(i,:))]);
write_text_file([output_dir '/matrix.' prefix '.E.N' num2str(Nkmer) '.txt'],[Sk(i) num2cell(Ek(i,:))]);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
Pc = corr(Pk');
Pe = corr(Ek');
[~,i] = cluster_sort(Pc+Pe);
imagesc(Pc(i,i),[-1 1]);
axis square;
colormap(gene_colormap(1));
N = Sk(i);
set(gca,'xtick',1:size(N,1),'xticklabel',N,'ytick',1:size(N,1),'yticklabel',N);
saveas(h, [output_dir '/matrix.' prefix '.N' num2str(Nkmer) '.corr.jpg'],'jpg');
