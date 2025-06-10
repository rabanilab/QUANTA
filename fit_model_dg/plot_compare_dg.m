function h = plot_compare_dg(x,y,zid,all_names,orgid,drange,orgid2)

if (nargin < 7)
    orgid2 = 'zfish';
end

pctl = 2;
fold = log2(1.5);

i1 = (x>=prctile(x,pctl));%.*(x<=prctile(x,100-pctl)) == 1;
i2 = (y>=prctile(y,pctl));%.*(y<=prctile(y,100-pctl)) == 1;
i = i1.*i2 == 1;

f = mean(y(i)-x(i));
p = prctile(abs(x(i)-(y(i)-f)),25);


h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(1,2,1);
hold on;
dscatter(x(i),y(i),'MSIZE',60);
line(drange,drange,'LineStyle',':','color','k','linewidth',1);
line(drange,drange+f,'LineStyle','-','color','k','linewidth',1);
line(drange,drange+f-p,'LineStyle','-','color','r','linewidth',1);
line(drange,drange+f+p,'LineStyle','-','color','r','linewidth',1);
if (~isempty(all_names))
    for j = 1:max(size(all_names))
        k = strcmp(zid,all_names{j});
        plot(x(k),y(k),'ok','markersize',10);
        text(x(k),y(k),all_names{j},'fontsize',18);
    end
end
hold off;
axis square;
axis tight;
set(gca,'xlim',drange,'ylim',drange,'xtick',-20:20,'ytick',-20:20,'fontsize',18);
xlabel(sprintf('%s (dg rate, log2)', orgid));
ylabel([orgid2 ' (dg rate, log2)']);
title(sprintf('dg rate, log2 (n=%d, r=%.2f, f=%.1f)',sum(i),corr(x(i),y(i)),2.^f));

subplot(1,2,2);
hold on;
dscatter(x(i),y(i)-f,'MSIZE',60);
line(drange,drange,'LineStyle',':','color','k','linewidth',1);
line(drange,drange-fold,'LineStyle','-','color','r','linewidth',1);
line(drange,drange+fold,'LineStyle','-','color','r','linewidth',1);
if (~isempty(all_names))
    for j = 1:max(size(all_names))
        k = strcmp(zid,all_names{j});
        plot(x(k),y(k)-f,'ok','markersize',10);
        text(x(k),y(k)-f,all_names{j},'fontsize',18);
    end
end
q = sum(abs(x(i)-(y(i)-f)) <= fold);
hold off;
axis square;
axis tight;
set(gca,'xlim',drange,'ylim',drange,'xtick',-20:20,'ytick',-20:20,'fontsize',18);
xlabel(sprintf('%s (dg rate, log2)', orgid));
ylabel([orgid2 ' (dg rate, log2) - scaled']);
title(sprintf('dg rate, log2 (n=%d, %d (%.1f%%))',sum(i),q,100*q./sum(i)));
