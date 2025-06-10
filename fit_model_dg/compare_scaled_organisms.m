function [mtid,mtidz,LLR,Q] = compare_scaled_organisms(ofile,opref,kmer_dir,max_scaleT,all_data_r,all_data_a,data_pref)

if (nargin < 7)
    data_pref = 'data_from_drive/workdir/exp_df_classified_';
end

kmer_range = 4:7;
kmer_alpha = 0.05;
kmer_esize = 30;

kmer_dir_z = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/gene_kmers';
all_data_r_z = {'Medina2021_ribo' 'Zhao2017' 'Meier2018'};% 'Vejnar2019'};
all_data_a_z = {};%{'Medina2021_polya' 'Yang2019' 'Pauli2011' 'Harvey2013'};

% ---------------------------------------------------------------------
% current organism expression data
% ---------------------------------------------------------------------
nr = max(size(all_data_r));
na = max(size(all_data_a));
load('maternal_param.mat','mtid');
nmtid = size(mtid,1);

mT = cell(1,nr);
mD = cell(1,nr);
for j = 1:nr
    all_data_r{j}
    [D,t,iM,iP,gid,tid,gc] = load_data([data_pref all_data_r{j} '.txt']);
    [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = D(k2(k2>0),iM);
    mD{j} = nan(nmtid,size(y,2));
    mD{j}(k1,:) = y;
    mT{j} = x;
end
aT = cell(1,na);
aD = cell(1,na);
for j = 1:na
    all_data_a{j}
    [D,t,iM,iP,gid,tid,gc] = load_data([data_pref all_data_a{j} '.txt']);
    [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = D(k2(k2>0),iM);
    aD{j} = nan(nmtid,size(y,2));
    aD{j}(k1,:) = y;
    aT{j} = x;
end

all_data = [all_data_r all_data_a];
D = [mD aD];
T = [mT aT];

% ---------------------------------------------------------------------
% zebrafish expression data
% ---------------------------------------------------------------------
nr = max(size(all_data_r_z));
na = max(size(all_data_a_z));
load('../1_zebrafish_Nov2023/maternal_param.mat','mtid');
mtidz = mtid;
nmtid = size(mtid,1);
load('maternal_param.mat','mtid');

mTz = cell(1,nr);
mDz = cell(1,nr);
for j = 1:nr
    all_data_r_z{j}
    [X,t,iM,iP,gid,tid,gc] = load_data(['../1_zebrafish_Nov2023/data_from_drive/workdir/exp_df_classified_' all_data_r_z{j} '.txt']);
    [k1,k2] = ismember(mtidz,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = X(k2(k2>0),iM);
    mDz{j} = nan(nmtid,size(y,2));
    mDz{j}(k1,:) = y;
    mTz{j} = x;
end
aTz = cell(1,na);
aDz = cell(1,na);
for j = 1:na
    all_data_a_z{j}
    [X,t,iM,iP,gid,tid,gc] = load_data(['../1_zebrafish_Nov2023/data_from_drive/workdir/exp_df_classified_' all_data_a_z{j} '.txt']);
    [k1,k2] = ismember(mtidz,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = X(k2(k2>0),iM);
    aDz{j} = nan(nmtid,size(y,2));
    aDz{j}(k1,:) = y;
    aTz{j} = x;
end

all_data_z = [all_data_r_z all_data_a_z];
Dz = [mDz aDz];
Tz = [mTz aTz];

% ---------------------------------------------------------------------
% intersect ids
% ---------------------------------------------------------------------
f = fopen(ofile);
X = textscan(f,'%s %s');
fclose(f);

[~,j1,j2] = intersect(mtid,X{2}); %[mtid(j1) X{2}(j2(j2>0))]
for j = 1:max(size(D))
    D{j} = D{j}(j1,:);
end
mtid = mtid(j1);
X{1} = X{1}(j2(j2>0));
X{2} = X{2}(j2(j2>0));

[~,jz1,jz2] = intersect(mtidz,X{1}); %[mtidz(jz1) X{1}(jz2(jz2>0))]
for j = 1:max(size(D))
    D{j} = D{j}(jz2(jz2>0),:);
end
mtid = mtid(jz2(jz2>0));
for j = 1:max(size(Dz))
    Dz{j} = Dz{j}(jz1,:);
end
mtidz = mtidz(jz1);

% ---------------------------------------------------------------------
% scale time
% scaled Ti = Ti * (maxT / max(Ti))
% ---------------------------------------------------------------------
mxT = 0;
for i = 1:max(size(T))
    mxT = max(mxT,max(T{i}));
end

for i = 1:max(size(T))
    all_data{i}
    T{i}    
    T{i} = T{i}.*(max_scaleT./mxT);
    T{i}
end

% ---------------------------------------------------------------------
% statistical test
% ---------------------------------------------------------------------
llr_alpha = 0.01;

Q = [];
ID = {};
mkdir([opref '.llr']);
for i = 1:max(size(all_data_z))
    for j = 1:max(size(all_data))
        ID = [ID;{[all_data{j} '.' all_data_z{i}]}];
        prefname = [opref '.llr/' all_data{j} '.' all_data_z{i}];
        prefname
        Qi = scaled_fit(prefname,mtidz,[T(j) Tz(i)],[D(j) Dz(i)],[all_data(j) all_data_z(i)]);
        Q = [Q Qi];
    end
end

% FDR correction to account for a different number of datasets in each
% organism
Qx = reshape(mafdr(Q(:),'BHFDR',true),size(Q,1),size(Q,2));
LLR = Qx < llr_alpha;
W = [['ensid' 'gid' ID'];[mtidz mtid num2cell(-1*log10(Qx))]];
write_text_file([opref '.llrp.txt'],W);

[u,~,t] = unique(LLR,'rows');
c = accumarray(t,1);
ui = regexprep(cellstr(num2str(u)),'  *',',');
W = sortrows([ui num2cell([c round(100*c./sum(c),1)])],2)
write_text_file([opref '.llr.txt'],W);

W = [{'id' 'reject' 'pct'};[ID num2cell([sum(LLR)' 100*sum(LLR)'/size(LLR,1)])]]
write_text_file([opref '.cnt.txt'],W);

% select scaled genes
k = sum(LLR,2)./size(LLR,2) > 1/3;
W = [{'reject' 'retain'}; num2cell([sum(k==1) sum(k==0)]); num2cell(round(100*[sum(k==1) sum(k==0)]./size(k,1),1))]
write_text_file([opref '.sum.txt'],W);
write_text_file([opref '.reject.txt'],[mtidz mtid num2cell([k sum(LLR,2) sum(LLR,2)./size(LLR,2)])]);

% ---------------------------------------------------------------------
% kmer enrichments
% ---------------------------------------------------------------------
mkdir([opref '.kmers']);

volcano_hyg(mtidz,k,[opref '.kmers/z'],kmer_range,kmer_esize,kmer_alpha,kmer_dir_z,[-100 100]);
volcano_hyg(mtid,k,[opref '.kmers/o'],kmer_range,kmer_esize,kmer_alpha,kmer_dir,[-100 100]);
% for i = 1:size(LLR,2)
%     volcano_hyg(mtidz,LLR(:,i),[opref '.kmers/z.' ID{i}],kmer_range,kmer_esize,kmer_alpha,kmer_dir_z,[-100 100]);
%     volcano_hyg(mtid,LLR(:,i),[opref '.kmers/o.' ID{i}],kmer_range,kmer_esize,kmer_alpha,kmer_dir,[-100 100]);
% end

if (exist('../scale_zfish.txt','file'))
    X = importdata('../scale_zfish.txt');
    volcano_hyg(X.textdata(:,1),X.data==0,[opref '.kmers/zall'],kmer_range,kmer_esize,kmer_alpha,kmer_dir_z,[-100 100]);
end


function Q = scaled_fit(result_pref,mtid,mT,mD,all_data_r)
% mP1 mR1 mE1 = single fit parameters

% single model fit
if (exist([result_pref '.model.1p.mat'],'file'))
    load([result_pref '.model.1p.mat'],'mtid','mP','mR','mE');
else
    [mP,mR,mE] = run_fit(mT,mD,[]);
    save([result_pref '.model.1p.mat'],'mtid','mP','mR','mE');
end

% scaled model fit
if (exist([result_pref '.model.tp.mat'],'file'))
    load([result_pref '.model.tp.mat'],'mtid','mPt','mRt','mEt');
else
    [mPt,mRt,mEt] = run_fit_temp(mT,mD,[]);
    save([result_pref '.model.tp.mat'],'mtid','mPt','mRt','mEt');
end

% likelihood ratio test
% LLR = 1: reject L1 (single model) for L2 (two models)
Q = model_llr_temp(mT,mD,mP,mPt);


% plot comparison results
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
LLR = Q < 0.01;
d1 = (log(2)./mP{1}(:,2) - log(2)./mP{2}(:,2)); % half-life
d2 = (mP{1}(:,3) - mP{2}(:,3)); % onset time
k = ones(size(mtid,1),1);
for i = 1:max(size(mR))
    k = k.*(mR{i}>=0.7).*(mE{i}<=20).*(mP{i}(:,1)>=3) == 1;
end

subplot(2,2,1);
y = hist(LLR,[0 1]);
b = bar(y./sum(y));
axis tight;
ylabel('fraction of genes');
set(gca,'xticklabel',{'L1' 'L2'},'ylim',[0 1],'fontsize',14);
b.FaceColor = 'flat';
for i = 1:2
    text(i-0.3,y(i)./sum(y)-0.05,sprintf('n=%d\n(%.1f%%)',y(i),100*y(i)./sum(y)),'fontsize',14);
    b.CData(i,:) = [0.8 0.8 0.8];
end
subplot(2,2,2);
y = hist(LLR(k),[0 1]);
b = bar(y./sum(y));
axis tight;
ylabel('fraction of genes');
set(gca,'xticklabel',{'L1' 'L2'},'ylim',[0 1],'fontsize',14);
b.FaceColor = 'flat';
for i = 1:2
    text(i-0.3,y(i)./sum(y)-0.05,sprintf('n=%d\n(%.1f%%)',y(i),100*y(i)./sum(y)),'fontsize',14);
    b.CData(i,:) = [0.8 0.8 0.8];
end
title('(R2>=0.8)');
subplot(2,2,3);
x = -7:0.2:7;
y = [];
y(1,:) = hist(d1(LLR==0),x);
y(2,:) = hist(d1(LLR==1),x);
plot(x,y'./sum(y'),'-','linewidth',2);
xlabel(sprintf('difference in half-life (%s-%s)',all_data_r{1},all_data_r{2}));
ylabel('fraction of genes');
axis tight;
legend({'L1' 'L2'},'Box','off');
set(gca,'fontsize',14);
subplot(2,2,4);
x = -7:0.2:7;
y = [];
y(1,:) = hist(d2(LLR==0),x);
y(2,:) = hist(d2(LLR==1),x);
plot(x,y'./sum(y'),'-','linewidth',2);
xlabel(sprintf('difference in onset (%s-%s)', ...
    regexprep(all_data_r{1},'_',' '),regexprep(all_data_r{2},'_',' ')));
ylabel('fraction of genes');
axis tight;
legend({'L1' 'L2'},'Box','off');
set(gca,'fontsize',14);
saveas(h, [result_pref '.llr.jpg'],'jpg');
saveas(h, [result_pref '.llr.eps'],'epsc');

% if (sum(k) > 1)
%     h = plot_compare(mP{1},mP{2},k,all_data_r(1:2));
%     saveas(h,[result_pref '.model_param.jpg'],'jpg');
%     h = plot_compare(mP{1},mPt{1},k,{all_data_r{1} 'joint'});
%     saveas(h,[result_pref '.model_param.' all_data_r{1} '.jpg'],'jpg');
%     h = plot_compare(mP{2},mPt{2},k,{all_data_r{2} 'joint'});
%     saveas(h,[result_pref '.model_param.' all_data_r{2} '.jpg'],'jpg');
% end
% y = (LLR==1).*k == 1;
% if (sum(y) > 1)
%     h = plot_compare(mP{1},mP{2},y,all_data_r(1:2));
%     saveas(h,[result_pref '.model_param.reject.jpg'],'jpg');
% end
% z = (LLR==0).*k == 1;
% if (sum(z) > 1)
%     h = plot_compare(mP{1},mP{2},z,all_data_r(1:2));
%     saveas(h,[result_pref '.model_param.retain.jpg'],'jpg');
% end
close all;



function [mP,mR,mE] = run_fit_temp(mT,mD,mzt)

n = max(size(mT));
m = size(mD{1},1);
fprintf('TMP model: fitting %d samples, %d reportrs\n',n,m);

% fit each reporter
P = nan(m,n+3);
R = nan(m,1);
E = nan(m,1);
parfor i = 1:m
    di = [];
    for j = 1:n
        di{j} = mD{j}(i,:);
    end
    [P(i,:),E(i,1),~,R(i,1)] = fit_model_temp(mT,di,[]);%mzt{1});
end

% final models
mP = cell(1,n);
mR = cell(1,n);
mE = cell(1,n);
for i = 1:n
    mP{i} = P(:,[i n+1:n+3]);
    mR{i} = R;
    mE{i} = E;
end


function [mP,mR,mE] = run_fit(mT,mD,mzt)

n = max(size(mT));
m = size(mD{1},1);
fprintf('SINGLE model: fitting %d samples, %d reportrs\n',n,m);

% fit each reporter
mP = cell(1,n);
mR = cell(1,n);
mE = cell(1,n);
parfor j = 1:n
    mP{j} = nan(m,4);
    mR{j} = nan(m,1);
    mE{j} = nan(m,1);
    for i = 1:m
        [mP{j}(i,:),mE{j}(i,1),~,mR{j}(i,1)] = fit_model(mT{j},mD{j}(i,:),[]);%mzt{1});
    end
end


function Q = model_llr_temp(mT,mD,mP,mPt)

sig = 0.5;

% load data
n = sum(cellfun(@isempty,mD)==0);

% estimate variance
sig2 = [];
for j = 1:n
    Y = mD{j};
    Y1 = Y(:,1:end-1);
    Y2 = Y(:,2:end);
    sig2(j) = sum(sum((0.5*(Y1-Y2)).^2))./(size(Y1,1)*size(Y1,2));
end
siga = sqrt(sig2) + sig;
fprintf('sigma = \n');
num2cell(siga)
% siga = 2*sqrt(sig2);
% siga = 2*repmat(sig,1,n);

% Model #1 [x0,dg,t0,t1]
for j = 1:n
    Y = mD{j};
    P1 = mPt{j};
    X1 = [];
    for i = 1:size(P1,1)
        X1(i,:) = dg_eval_model(mT{j},P1(i,:));
    end
    L1(:,j) = sum(log2(normpdf(Y,X1,siga(j))),2);
end
L1 = sum(L1,2);

% Model #2 [x0,da,dg,t0,t1]
for j = 1:n
    Y = mD{j};
    P2 = mP{j};
    X2 = [];
    for i = 1:size(P2,1)
        X2(i,:) = dg_eval_model(mT{j}, P2(i,:));
    end
    L2(:,j) = sum(log2(normpdf(Y,X2,siga(j))),2);
end
L2 = sum(L2,2);

fprintf('rows with L1>L2: %d\n', sum(L2<L1));
[~,Q] = lratiotest(L2,min(L1,L2),3);%n*size(P2,2)-size(P1,2));


function h = plot_compare(mP1,mP2,k,ids)

x = [mP1(k,1) mP2(k,1)];
d = [mP1(k,2) mP2(k,2)];
l = log(2)./[mP1(k,2) mP2(k,2)];
t = [mP1(k,3) mP2(k,3)];

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(2,2,1);
dscatter(d(:,1),d(:,2),'MSIZE',25);
set(gca,'xlim',[0 3],'ylim',[0 3],'fontsize',18);
axis square;
xlabel(ids{1});
ylabel(ids{2});
r2 = median(d(:,1))./median(d(:,2));
title(sprintf('degradation rate (ratio=%.2f, n=%d)',r2,size(d,1)));
hold on;
p2 = polyfit(d(:,1),d(:,2),1);
q = polyval(p2,[0 3]);
plot([0 3],[0 3],'-k');
plot([0 3],q,'-r');
hold off;

subplot(2,2,2);
dscatter(l(:,1),l(:,2),'MSIZE',25);
set(gca,'xlim',[0 5],'ylim',[0 5],'fontsize',18);
axis square;
xlabel(ids{1});
ylabel(ids{2});
r2 = median(l(:,1))./median(l(:,2));
title(sprintf('half-life (ratio=%.2f, n=%d)',r2,size(l,1)));
hold on;
p2 = polyfit(l(:,1),l(:,2),1);
q = polyval(p2,[0 5]);
plot([0 5],[0 5],'-k');
plot([0 5],q,'-r');
hold off;

subplot(2,2,3);
dscatter(x(:,1),x(:,2),'MSIZE',25);
set(gca,'xlim',[4 12],'ylim',[4 12],'fontsize',18);
axis square;
xlabel(ids{1});
ylabel(ids{2});
r1 = median(x(:,1)) - median(x(:,2));
title(sprintf('logX0 (ratio=%.2f)',2.^r1));
hold on;
p1 = polyfit(x(:,1),x(:,2),1);
q = polyval(p1,[4 12]);
plot([4 12],[4 12],'-k');
plot([4 12],q,'-r');
hold off;

subplot(2,2,4);
dscatter(t(:,1),t(:,2),'MSIZE',25);
set(gca,'xlim',[0 11],'ylim',[0 11],'fontsize',18);
axis square;
xlabel(ids{1});
ylabel(ids{2});
r3 = median(t(:,1))./median(t(:,2));
title(sprintf('onset time (ratio=%.2f)',r3));
hold on;
m = median(t(:,1));
line([m m],[0 11],'LineStyle','-','color','k','linewidth',1);
m = median(t(:,2));
line([0 11],[m m],'LineStyle','-','color','k','linewidth',1);
