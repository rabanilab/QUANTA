function h = plot_compare_rates(nx,orgid,f,drange,suffix,orgid2)

if (nargin < 5)
    suffix = '';
end
if (nargin < 6)
    orgid2 = 'zfish';
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(2,4,1);
X = importdata(['orthologies.dg' suffix '.txt']);
i = ismember(X.textdata(:,1),nx);
x = log2(X.data(i,2));
y = log2(X.data(i,1));
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);
dscatter(x,y,'MSIZE',20);
hold on;
line(drange,drange,'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
axis tight;
ylabel(sprintf('%s (dg rate, log2)',orgid));
xlabel([orgid2 ' (dg rate, log2)']);
[c,p] = corr(x,y);
title(sprintf('dg (n = %d, r = %.2f, p < %.0e)',size(x,1),c,p));
set(gca,'fontsize',14);

subplot(2,4,5);
X = importdata(['orthologies.dg' suffix '.txt']);
i = ismember(X.textdata(:,1),nx);
x = X.data(i,2);
%x = (x - min(x))./(max(x) - min(x));
y = X.data(i,1);
y = y./f;
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);
dscatter(x,y,'MSIZE',20);
hold on;
mn = min([x;y]);
mx = max([x;y]);
line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
axis tight;
ylabel(sprintf('%s (dg rate) - scaled',orgid));
xlabel([orgid2 ' (dg rate)']);
[c,p] = corr(x,y);
title(sprintf('dg (n = %d, r = %.2f, p < %.0e)',size(x,1),c,p));
set(gca,'fontsize',14);

subplot(2,4,6);
X = importdata(['orthologies.t0' suffix '.txt']);
i = ismember(X.textdata(:,1),nx);
x = X.data(i,2);
%x = (x - min(x))./(max(x) - min(x));
y = X.data(i,1);
y = y*f;
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);
dscatter(x,y,'MSIZE',20);
hold on;
mn = min([x;y]);
mx = max([x;y]);
line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
axis tight;
ylabel(sprintf('%s (t0) - scaled',orgid));
xlabel([orgid2 ' (t0)']);
[c,p] = corr(x,y);
title(sprintf('t0 (n = %d, r = %.2f, p < %.0e)',size(x,1),c,p));
set(gca,'fontsize',14);

subplot(2,4,7);
X = importdata(['orthologies.da' suffix '.txt']);
i = ismember(X.textdata(:,1),nx);
x = (X.data(i,2));
%x = (x - min(x))./(max(x) - min(x));
y = (X.data(i,1));
y = y./f;
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);
dscatter(x,y,'MSIZE',20);
hold on;
mn = min([x;y]);
mx = max([x;y]);
line([mn mx],[mn mx],'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
axis tight;
ylabel(sprintf('%s (da rate) - scaled',orgid));
xlabel([orgid2 ' (da rate)']);
[c,p] = corr(x,y);
title(sprintf('da (n = %d, r = %.2f, p < %.0e)',size(x,1),c,p));
set(gca,'fontsize',14);

subplot(2,4,8);
X = importdata(['orthologies.x0' suffix '.txt']);
i = ismember(X.textdata(:,1),nx);
x = X.data(i,2);
%x = (x - min(x))./(max(x) - min(x));
y = X.data(i,1);
%y = (y - min(y))./(max(y) - min(y));
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);
dscatter(x,y,'MSIZE',20);
hold on;
line([-2 12],[-2 12],'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
axis tight;
ylabel(sprintf('%s (x0, log2)',orgid));
xlabel([orgid2 ' (x0, log2)']);
[c,p] = corr(x,y);
title(sprintf('x0 (n = %d, r = %.2f, p < %.0e)',size(x,1),c,p));
set(gca,'fontsize',14);
