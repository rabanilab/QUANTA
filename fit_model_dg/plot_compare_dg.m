function h = plot_compare_dg(x,y,zid,all_names,orgid,drange,orgid2)

if (nargin < 7)
    orgid2 = 'zfish';
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

f = mean(y-x);
p = prctile(abs(x-(y-f)),25);

hold on;
dscatter(x,y,'MSIZE',60);
line(drange,drange,'LineStyle',':','color','k','linewidth',1);
line(drange,drange+f,'LineStyle','-','color','k','linewidth',1);
line(drange,drange+f-p,'LineStyle','-','color','r','linewidth',1);
line(drange,drange+f+p,'LineStyle','-','color','r','linewidth',1);
if (~isempty(all_names))
    for j = 1:max(size(all_names))
        k = strcmp(zid,all_names{j});
        plot(x(k),y(k),'ok','markersize',10);
        text(x(k),y(k),all_names{j},'fontsize',35);
    end
end
hold off;
axis square;
axis tight;
set(gca,'xlim',drange,'ylim',drange,'xtick',-20:20,'ytick',-20:20,'fontsize',35);
xlabel(sprintf('%s (dg rate, log2)', orgid));
ylabel([orgid2 ' (dg rate, log2)']);
title(sprintf('dg rate, log2 (n=%d, r=%.2f, f=%.1f)',size(x,1),corr(x,y),2.^f));
