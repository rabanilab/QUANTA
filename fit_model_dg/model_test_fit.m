function [Xr2,Xe] = model_test_fit(all_data,mtid,mP,mE,mR,model,output_pref,minR2,maxE,data_id,pref)

if (nargin < 8)
    minR2 = 0.7;
end
if (nargin < 9)
    maxE = 7;
end
if (nargin < 10)
    data_id = [];
end
if (nargin < 11)
    pref = 'data_from_drive/exp_df_classified_';
end

nd = size(all_data,2);

% r-squared
fprintf('R2 > %.1f\n',minR2);
R2 = []; 
for j = 1:nd
    R2(:,j) = mR{j};
end
Xr2 = [[sum(R2>minR2,1)' sum(R2>minR2,1)'./size(R2,1)];[sum(R2(:)>minR2,1) sum(R2(:)>minR2,1)./size(R2(:),1)]];
[[all_data';'all'] num2cell(Xr2)]

% mean square error
fprintf('MSE < %.1f\n',maxE);
MSE = []; 
for j = 1:nd
    MSE(:,j) = mE{j};
end
Xe = [[sum(MSE<maxE,1)' sum(MSE<maxE,1)'./size(MSE,1)];[sum(MSE(:)<maxE,1) sum(MSE(:)<maxE,1)./size(MSE(:),1)]];
[[all_data';'all'] num2cell(Xe)]


h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
x = -1:0.05:1;
y = hist(R2,x);
plot(x,y./sum(y),'-','linewidth',2);
box off;
axis tight;
%set(gca,'ylim',[0 0.35]);
legend(regexprep(all_data,'_',' '),'box','off','location','bestOutside');
xlabel('r-squared');
ylabel('fraction of genes');
saveas(h, [output_pref '.rsq_genes.jpg'],'jpg');

clf;
x = -10:0.5:10;
y = hist(log2(MSE),x);
plot(x,y./sum(y),'-','linewidth',2);
box off;
axis tight;
%set(gca,'ylim',[0 0.2]);
legend(regexprep(all_data,'_',' '),'box','off','location','bestOutside');
xlabel('mean squared error; log2');
ylabel('fraction of genes');
saveas(h, [output_pref '.mse_genes.jpg'],'jpg');


% correlation data to predictions
[~,mD,mM,mM0] = dg_eval_all(all_data,mtid,mP,model,data_id,pref);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
xlim = [-4 12];
nk = ceil(nd/2);
r = [];
for j = 1:nd
    subplot(2,nk,j);
    x = mD{j}(:);
    x(x<xlim(1)) = xlim(1);
    x(x>xlim(2)) = xlim(2);
    y = mM{j}(:);
    y(y<xlim(1)) = xlim(1);
    y(y>xlim(2)) = xlim(2);
    i = (~isnan(x)).*(~isnan(y)) == 1;
    nx = 0.01*randn(sum(i),1);
    ny = 0.01*randn(sum(i),1);
    dscatter(x(i)+nx,y(i)+ny);
    axis square;
    set(gca,'xlim',xlim,'ylim',xlim,'fontsize',14,'box','on');
    xlabel('data');
    ylabel('model');
    r(j) = corr(x(i),y(i));
    title(sprintf('%s (r2=%.2f)',regexprep(all_data{j},'_',' '),r(j).^2));
end
saveas(h, [output_pref '.corr_datasets.jpg'],'jpg');

clf;
bar(r);
axis tight;
set(gca,'xticklabel',regexprep(all_data,'_',' '),'ylim',[0 1.2],'fontsize',16);
for j = 1:nd
    text(j,r(j),sprintf('%.2f',r(j)),'FontSize',18);
end
ylabel('r-squared model to data');
saveas(h, [output_pref '.corr_datasets.bar.jpg'],'jpg');
saveas(h, [output_pref '.corr_datasets.bar.eps'],'epsc');

clf;
x = [];
y = [];
z = [];
for j = 1:nd
    x = [x;mD{j}(:)];
    y = [y;mM{j}(:)];
    z = [z;mM0{j}(:)];
end
subplot(1,2,1);
i = (~isnan(x)).*(~isnan(y)) == 1;
dscatter(x(i),y(i));
axis square;
set(gca,'xlim',xlim,'ylim',xlim,'fontsize',14,'box','on');
xlabel('data');
ylabel('model');
r = corr(x(i),y(i));
title(sprintf('all data (r2=%.2f)',r.^2));
subplot(1,2,2);
i = (~isnan(x)).*(~isnan(z)) == 1;
dscatter(x(i),z(i));
axis square;
set(gca,'xlim',xlim,'ylim',xlim,'fontsize',14,'box','on');
xlabel('data');
ylabel('null model');
r = corr(x(i),z(i));
title(sprintf('all data (r2=%.2f)',r.^2));
saveas(h, [output_pref '.corr_all.jpg'],'jpg');

close all;

% compare between predictions on different datasets
plot_corr_all(all_data,mP,[output_pref '.corr_param']);




function plot_corr_all(all_data,P,opref)

i1 = [];
i2 = [];
for i = 1:max(size(all_data))
    for j = i+1:max(size(all_data))
        i1 = [i1 i];
        i2 = [i2 j];
    end
end
maxN = max(ceil(max(size(i1))/3),1);

nParam = size(P{1},2);

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
k = 1;

% initial levels
clf;
for j = 1:max(size(i1))
    x1 = P{i1(j)}(:,k);
    x2 = P{i2(j)}(:,k);
    plot_corr(x1,all_data{i1(j)},x2,all_data{i2(j)},j,maxN);
end
saveas(h, [opref '.x0.jpg'],'jpg');
k = k+1;

% deadenylation rate
if (nParam == 5)
    clf;
    for j = 1:max(size(i1))
        x1 = P{i1(j)}(:,k);
        x1(x1==0) = NaN;
        x2 = P{i2(j)}(:,k);
        x2(x2==0) = NaN;
        plot_corr(x1,all_data{i1(j)},x2,all_data{i2(j)},j,maxN);
    end
    saveas(h, [opref '.da.jpg'],'jpg');
    k = k+1;
end

% degradation rate
clf;
for j = 1:max(size(i1))
    x1 = log2(P{i1(j)}(:,k));
    x2 = log2(P{i2(j)}(:,k));
    x1(x1<log2(log(2)/8))=log2(log(2)/8);
    x2(x2<log2(log(2)/8))=log2(log(2)/8);
    plot_corr(x1,all_data{i1(j)},x2,all_data{i2(j)},j,maxN);
end
saveas(h, [opref '.dg.jpg'],'jpg');

% half-life
clf;
for j = 1:max(size(i1))
    x1 = log(2)./P{i1(j)}(:,k);
    x2 = log(2)./P{i2(j)}(:,k);
    plot_corr(x1,all_data{i1(j)},x2,all_data{i2(j)},j,maxN);
end
saveas(h, [opref '.hl.jpg'],'jpg');
k = k+1;

% onset time
clf;
for j = 1:max(size(i1))
    x1 = P{i1(j)}(:,k);
    x2 = P{i2(j)}(:,k);
    plot_corr(x1,all_data{i1(j)},x2,all_data{i2(j)},j,maxN);
end
saveas(h, [opref '.t0.jpg'],'jpg');



function plot_corr(x1,t1,x2,t2,j,maxN)

if (nargin < 6)
    maxN = 3;
end

subplot(3,maxN,j);
k = (~isnan(x1)).*(~isnan(x2)).*(~isinf(x1)).*(~isinf(x2))==1;
if (sum(k) > 10)
    rn1 = 0.01*randn(size(x1(k)));
    rn2 = 0.01*randn(size(x2(k)));
    dscatter(x1(k)+rn1,x2(k)+rn2,'MSIZE',20);
    minX = min([x1(k);x2(k)]);
    maxX = max([x1(k);x2(k)]);
    if (maxX-minX<0.1)
        maxX = maxX+1;
    end
    set(gca,'xlim',[minX maxX],'ylim',[minX maxX]);
    r = unique(round(minX:maxX));
    if (max(size(r))>5)
        dr = round(max(size(r))/5);
        if (dr<1)
            dr = 1;
        end
        r = r(1:dr:end);
    end
    set(gca,'xtick',r,'ytick',r,'fontsize',14);
    title(sprintf('r=%.2f, n=%d',corr(x1(k),x2(k)),sum(k)));
end
xlabel(regexprep(t1,'_',' '));
ylabel(regexprep(t2,'_',' '));
axis square;

