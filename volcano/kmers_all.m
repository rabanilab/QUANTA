function kmers_all(Dids,Dseq,Dval,Dnames,output_dir,Krange,Kdir,esize,palpha,erange,prange,sort_type,norm_len)

if (nargin < 6)
    Krange = 3:7;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    esize = 5;
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
if (nargin < 12)
    sort_type = 1; % 1 = by pvalues; 2 = by kmers
end
if (nargin < 13)
    norm_len = 1;
end

n = size(Dval,2);
minP = 1e-200;
min_overlap = min(Krange);
mxP = prange(2);
mxE = max(abs(erange));
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
if (exist([output_dir '/ks.mat'],'file'))
    load([output_dir '/ks.mat'],'W1','W2','Ethr');
    fprintf('load KS-test data (all significant kmers)\n');
    [Dnames' W1' W2' num2cell(Ethr')]
else
    W1 = cell(1,n);
    W2 = cell(1,n);
    Ethr = zeros(1,n);
    for i = 1:n
        [num2cell(i) Dnames(i)]
        
        % remove extreme lengths: very long/short sequences skew the statistical test
        I = (~isnan(Dval(:,i))).*(Dlen > (1/len_fold)*median(Dlen)).*(Dlen < len_fold*median(Dlen)) == 1;

        % test if there is a correlation between length and property
        p = fit_linear_model(Dlen(I),Dval(I,i));
        x = [1 max(Dlen(I))];
        y = polyval(p,x);
        zy = y - (polyval(p,x)-p(2));
        zDval = Dval(I,i) - (polyval(p,Dlen(I)) - p(2));

        if (sum(I)>1)
            h = figure;
            subplot(1,2,1);
            hold on;
            w1 = 0.01*randn(size(Dlen(I)));
            w2 = 0.01*randn(size(Dlen(I)));
            dscatter(Dlen(I)+w1,Dval(I,i)+w2,'MSIZE',25);
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
            dscatter(Dlen(I)+w1,zDval+w2,'MSIZE',25);
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
                [W1{i},W2{i},Ethr(i)] = volcano_ks(Dids(I),zDval,[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
            else
                [W1{i},W2{i},Ethr(i)] = volcano_ks(Dids(I),Dval(I,i),[output_dir '/' Dnames{i}],Krange,esize,palpha,Kdir,erange,prange);
            end
        end
    end
    save([output_dir '/ks.mat'],'W1','W2','Ethr');
end

% write tables of significant kmers
S_all = [];
for i = 1:n
    if (~isempty(W1{i}))
        I = abs(cell2mat(W1{i}(:,5)))>=Ethr(i);%esize;
        S_all = [S_all;W1{i}(I,1)];
    end
    if (~isempty(W2{i}))
        I = abs(cell2mat(W2{i}(:,5)))>=Ethr(i);%esize;
        S_all = [S_all;W2{i}(I,1)];
    end
end
S_all = unique(S_all);
Q_all = ones(size(S_all,1),n);
E_all = zeros(size(S_all,1),n);
for i = 1:n
    if (~isempty(W1{i}))
        I = abs(cell2mat(W1{i}(:,5)))>=Ethr(i);%esize;
        wi = W1{i}(I,:);
        [w1,w2] = ismember(S_all,wi(:,1));
        if (sum(w1)>0)
            Q_all(w1,i) = cell2mat(wi(w2(w2>0),2));
            E_all(w1,i) = cell2mat(wi(w2(w2>0),5));
        end
    end
    if (~isempty(W2{i}))
        I = abs(cell2mat(W2{i}(:,5)))>=Ethr(i);%esize;
        wi = W2{i}(I,:);
        [w1,w2] = ismember(S_all,wi(:,1));
        if (sum(w1)>0)
            Q_all(w1,i) = cell2mat(wi(w2(w2>0),3));
            E_all(w1,i) = cell2mat(wi(w2(w2>0),5));
        end
    end
end
if (~isempty(S_all))
    P_all = -1*sign(E_all).*log10(Q_all);
else
    P_all = [];
end
[~,k] = sort(nanmax(abs(P_all),[],2),'descend');
write_text_file([output_dir '/ks.P.txt'],[['kmer' Dnames];[S_all(k) num2cell(P_all(k,:))]]);
write_text_file([output_dir '/ks.E.txt'],[['kmer' Dnames];[S_all(k) num2cell([E_all(k,:)])]]);
 

% ---------------------------------------------------------------------
% kmer peaks & motif prediction
% ---------------------------------------------------------------------
if (exist([output_dir '/pwm.mat'],'file'))
    load([output_dir '/pwm.mat'],'Q1','Q2','PWM1','PWM2');
    fprintf('load PWM data (all significant kmers)\n');
    [Q1;Q2]
else
    % negative weights
    Q1 = cell(n,1);
    PWM1 = cell(n,1);
    parfor i = 1:n
        fprintf('**** %s NEG ****\n',Dnames{i})
        Kid = W1{i}(:,1);
        Kweight = -1*log10(cell2mat(W1{i}(:,2)));
        I = abs(cell2mat(W1{i}(:,5)))>=Ethr(i);%esize;
        if (sum(I)==0)
            qi = [];
            pwmi = [];
        else
            [~,~,Z] = fit_peaks(Dids,Dseq,Kid(I),Kweight(I),[output_dir '/peaks.neg.' num2str(i)]);
            [qi,pwmi] = fit_pwms(Z);
            %Z = sortrows([num2cell(10.^Kweight(I)) Kid(I)],-1);
            %[qi,pwmi] = fit_pwms(Z);
        end
        Q1{i} = qi;
        PWM1{i} = pwmi;
    end

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    nc = 6;
    for i = 1:n
        pwmi = PWM1{i};
        nc = max(nc,max(size(pwmi)));
    end
    for i = 1:n
        qi = Q1{i};
        pwmi = PWM1{i};
        plot_pwm(n,nc,i,pwmi,qi,Dnames{i});
    end
    saveas(h,[output_dir '/pwm.neg.svg'],'svg');
    close all;

    % positive weights
    Q2 = cell(n,1);
    PWM2 = cell(n,1);
    parfor i = 1:n
        fprintf('**** %s POS ****\n',Dnames{i})
        Kid = W2{i}(:,1);
        Kweight = -1*log10(cell2mat(W2{i}(:,3)));
        I = abs(cell2mat(W2{i}(:,5)))>=Ethr(i);%esize;
        if (sum(I)==0)
            qi = [];
            pwmi = [];
        else
            [~,~,Z] = fit_peaks(Dids,Dseq,Kid(I),Kweight(I),[output_dir '/peaks.pos.' num2str(i)]);
            [qi,pwmi] = fit_pwms(Z);
            %Z = sortrows([num2cell(10.^Kweight(I)) Kid(I)],-1);
            %[qi,pwmi] = fit_pwms(Z);
        end
        Q2{i} = qi;
        PWM2{i} = pwmi;
    end
    save([output_dir '/pwm.mat'],'Q1','Q2','PWM1','PWM2');

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    nc = 6;
    for i = 1:n
        pwmi = PWM2{i};
        nc = max(nc,max(size(pwmi)));
    end
    for i = 1:n
        qi = Q2{i};
        pwmi = PWM2{i};
        plot_pwm(n,nc,i,pwmi,qi,Dnames{i});
    end
    saveas(h,[output_dir '/pwm.pos.svg'],'svg');
    close all;
end

% plot pwms
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

% negative peaks
clf;
nc = 6;
for i = 1:n
    pwmi = PWM1{i};
    nc = max(nc,max(size(pwmi)));
end
for i = 1:n
    qi = Q1{i};
    pwmi = PWM1{i};
    plot_pwm(n,nc,i,pwmi,qi,Dnames{i});
end
saveas(h,[output_dir '/pwm.neg.svg'],'svg');

% positive peaks
clf;
nc = 6;
for i = 1:n
    pwmi = PWM2{i};
    nc = max(nc,max(size(pwmi)));
end
for i = 1:n
    qi = Q2{i};
    pwmi = PWM2{i};
    plot_pwm(n,nc,i,pwmi,qi,Dnames{i});
end
saveas(h,[output_dir '/pwm.pos.svg'],'svg');

close all;


% ---------------------------------------------------------------------
% filter kmers: 
%  1. kmers that are contained within another kmer with a more significant p-value
%  2. kmers that contain another kmer with a more significant p-value
% P1/2 = [kmer] [pvalue] [effect size] [no. of kmers merged]
% ---------------------------------------------------------------------
if (exist([output_dir '/kmer.mat'],'file'))
    load([output_dir '/kmer.mat'],'P1','P2');
    fprintf('load filtered KMER data\n');
    [Dnames' P1' P2']
else
    P1 = cell(1,n);
    P2 = cell(1,n);
    for i = 1:n
        if (~isempty(W1{i}))
            score = cell2mat(W1{i}(:,[2 5]));
            I = abs(cell2mat(W1{i}(:,5)))>=Ethr(i);%esize;
            if (sum(I) > 0)
                P1{i} = filter_kmers_by_seq(W1{i}(I,1),score(I,:),min_overlap);
            end
        end
        if (~isempty(W2{i}))
            score = cell2mat(W2{i}(:,[3 5]));
            score(:,2) = -1*score(:,2);
            I = abs(cell2mat(W2{i}(:,5)))>=Ethr(i);%esize;
            if (sum(I) > 0)
                P2{i} = filter_kmers_by_seq(W2{i}(I,1),score(I,:),min_overlap);
                P2{i}(:,3) = num2cell(-1*cell2mat(P2{i}(:,3)));
            end
        end
    end
    save([output_dir '/kmer.mat'],'P1','P2');
end

% volcano plots
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

nc = 4;
nr = ceil(n/nc);
for j = 1:n
    P = [];
    E = [];
    S = [];
    if (~isempty(P1{j}))
        P = -1*log10(cell2mat(P1{j}(:,2)));
        P(P>300) = 300;% maximal pvalue is 1e-300; above that is 0 which is inf
        E = cell2mat(P1{j}(:,3));
        S = P1{j}(:,1);
    end
    if (~isempty(P2{j}))
        P = [P;-1*log10(cell2mat(P2{j}(:,2)))];
        P(P>300) = 300;% maximal pvalue is 1e-300; above that is 0 which is inf
        E = [E;cell2mat(P2{j}(:,3))];
        S = [S;P2{j}(:,1)];
    end

    subplot(nr,nc,j);
    if (~isempty(P))
        [k0,k1] = plot_volcano(P,E,S,Dnames{j},palpha,Ethr(j),0);%esize,0);
        if (sum(k0)>0)
            W = sortrows([S(k0) num2cell([10.^(-1*P(k0)) E(k0)])],[2 3]);
            write_text_file([output_dir '/kfiltered.' Dnames{j} '.volcano.txt'],W);
        end
        if (sum(k1)>0)
            W = sortrows([S(k1) num2cell([10.^(-1*P(k1)) E(k1)])],[2 3]);
            write_text_file([output_dir '/kfiltered.' Dnames{j} '.volcano.e.txt'],W);
        end
    else
        title(Dnames{j});
        set(gca,'fontsize',12);
        axis square;
    end
end
saveas(h, [output_dir '/kfiltered.volcano.jpg'],'jpg');
saveas(h, [output_dir '/kfiltered.volcano.svg'],'svg');
close all;


% ---------------------------------------------------------------------
% combine kmers from all different datasets
% ** Use only the identity of filtered k-mers, not p-values 
%    (take pvalues from the original set)
% ---------------------------------------------------------------------
if (exist([output_dir '/kmer_all.mat'],'file'))
    load([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all');
    fprintf('load combined KMER data\n');
    fprintf('total: %d kmers (all datasets)\n', size(S_all,1));
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
        P_all(P_all>300) = 300;% maximal pvalue is 1e-300; above that is 0 which is inf
    else
        P_all = [];
    end
    save([output_dir '/kmer_all.mat'],'S_all','E_all','Q_all','P_all');
end

% plot combined results
if (~isempty(S_all))

    % volcano plot
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    [~,mxJ] = max(abs(P_all),[],2);
    mP = zeros(size(P_all,1),1);
    mE = zeros(size(E_all,1),1);
    for i = 1:size(mxJ,1)
        mP(i,1) = P_all(i,mxJ(i));
        mE(i,1) = E_all(i,mxJ(i));
    end
    plot_volcano(mP,mE,S_all,'max',palpha,min(Ethr),1);%esize,1);
    saveas(h, [output_dir '/matrix.all.volcano.svg'],'svg');

    close all;

    % matrices
    if (~isempty(Dseq))
        P = P_all;
        P(E_all<0) = -1*P_all(E_all<0);
        P(isnan(P_all)) = 0;
        E = E_all;
        E(isnan(E_all)) = 0;

        % all kmers that pass both pvalue and esize thresholds
        k = zeros(size(P,1),1);
        for j = 1:n
            pj = (abs(P(:,j))>=-1*log10(palpha)).*(abs(E(:,j))>=Ethr(j)) == 1;%esize) == 1;
            k = k + pj > 0;
        end
        Sk = S_all(k,:);
        Ek = E(k,:);
        Pk = P(k,:);

        Nk = zeros(size(Sk,1),2);
        for i = 1:size(Sk,1)
            Nk(i,1) = sum(cellfun(@isempty,regexp(Dseq,Sk{i})) == 0);
            Nk(i,2) = sum(cellfun(@isempty,regexp(Dseq,Sk{i})) == 1);
        end

        if (sum(k) <= 500)
            %[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),max(Ethr),mxP,mxE);%esize,mxP,mxE);
            [h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,-1*log10(palpha),max(Ethr),mxP,mxE,[],1);%esize,mxP,mxE);
            set(gca,'xtick',[]);
            saveas(h, [output_dir '/matrix.all.jpg'],'jpg');
            write_text_file([output_dir '/matrix.P.all.txt'],[['kmer' Dnames 'withK' 'noK'];[Sk(i) num2cell([Pk(i,:) Nk(i,:)])]]);
            write_text_file([output_dir '/matrix.E.all.txt'],[['kmer' Dnames 'withK' 'noK'];[Sk(i) num2cell([Ek(i,:) Nk(i,:)])]]);
        else
            write_text_file([output_dir '/matrix.P.all.txt'],[['kmer' Dnames 'withK' 'noK'];[Sk num2cell([Pk Nk])]]);
            write_text_file([output_dir '/matrix.E.all.txt'],[['kmer' Dnames 'withK' 'noK'];[Sk num2cell([Ek Nk])]]);
        end

        % top N kmers per column
        plot_Nkmers(20,P,E,S_all,Dnames,sort_type,palpha,max(Ethr),mxP,mxE,output_dir);%esize,mxP,mxE,output_dir);
        plot_Nkmers(10,P,E,S_all,Dnames,sort_type,palpha,max(Ethr),mxP,mxE,output_dir);%esize,mxP,mxE,output_dir);

        % top 50 kmers per column (separate)
        for j = 1:n
            plot_Nkmers_col(j,50,P,E,S_all,Dnames,sort_type,palpha,max(Ethr),mxP,mxE,output_dir);%esize,mxP,mxE,output_dir);
        end
        close all;
    end
end


function p = fit_linear_model(x,y)

u1 = unique(x);
if (max(size(u1)) <= 2)
    p = [0 mean(y)];
else
    %p = polyfit(x,y,1);
    p = robustfit(x,y);
    p = p(end:-1:1)';
end

function plot_pwm(n,nc,i,pwmi,qi,iname)

for j = 1:max(size(pwmi))
    x = subplot(n,nc,j+(i-1)*nc);
    pwm_logo(pwmi{j},1000,'ACGU',x);
    xlabel('');
    ylabel('');
    set(gca,'xtick',[],'ytick',[],'fontsize',10);
    title(sprintf('%.0f%% (n=%d)',100*qi{j,3},qi{j,2}));
    if (j==1)
        ylabel(iname);
    end
end
if (max(size(pwmi))==0)
    subplot(n,nc,1+(i-1)*nc);
    box on;
    xlabel('');
    ylabel('');
    set(gca,'xtick',[],'ytick',[],'fontsize',10);
    ylabel(iname);
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

%[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE);
[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE,[],1);
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
%[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE,I);
[h,i] = plot_kmer_hitmap(Sk,Ek,Pk,Dnames,sort_type,[],esize,mxP,mxE,I,1);
set(gca,'xtick',[]);
saveas(h, [output_dir '/matrix.' Dnames{j} '.jpg'],'jpg');
write_text_file([output_dir '/matrix.P.' Dnames{j} '.txt'],[Sk(i) num2cell(Pk(i,:))]);
write_text_file([output_dir '/matrix.E.' Dnames{j} '.txt'],[Sk(i) num2cell(Ek(i,:))]);


function [k0,k1] = plot_volcano(P,E,S,name,palpha,esize,add_text)

k0 = isnan(P) + isnan(E) == 0;
k1 = (P >= -1*log10(palpha)).*(abs(E) >= esize) == 1;
mxPj = max(max(P),-1*log10(palpha)) + 1;
mxEj = max(max(abs(E)),esize) + 0.1;

hold on;

% % plot all points (grayscale)
% s = plot(E,P,'.k','markerSize',20);
% alpha(s,.5);

% write text
if (add_text)
    for i = find(k1)'
        seq = regexprep(S(i),'T','U');
        text(E(i),P(i),seq,'fontsize',10);
    end
end

% plot significant points
if (sum(k0) > 10)
    n1 = zeros(size(E(k0)));
    n2 = zeros(size(P(k0)));
    if (max(size(unique(E(k0))))<2)
        n1 = 0.01*randn(size(E(k0)));
    end
    if (max(size(unique(P(k0))))<2)
        n2 = 0.01*randn(size(P(k0)));
    end
    dscatter(E(k0)+n1,P(k0)+n2,'MSIZE',50,'MARKER','o');
else
    plot(E(k0),P(k0),'.','markerSize',20);
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
title(sprintf('%s (n=%d,n=%d)',name,sum(k0),sum(k1)));
set(gca,'xlim',[-1*mxEj mxEj],'ylim',[0 mxPj],'fontsize',12);
