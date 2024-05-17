function model_compare(P1,P2,Pids,prefix)
% P1 = degradation model fits
% P2 = polyA model fits

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

oid = regexprep(Pids,';','');
oid = regexprep(oid,' ','.');

n = max(size(P1));
for i = 1:n
    clf;
    subplot(2,2,1);
    plot_scatter(P1{i}(:,1),P2{i}(:,1),'initial RNA (log2)');
    subplot(2,2,2);
    plot_scatter(P1{i}(:,2),P2{i}(:,2),'dg vs. da rate (1/hr)');
    subplot(2,2,3);
    plot_scatter(P1{i}(:,2),P2{i}(:,3),'dg rate (1/hr)');
    subplot(2,2,4);
    plot_scatter(P1{i}(:,3),P2{i}(:,4),'t0 (hr)');
    saveas(h,[prefix '.' oid{i} '.jpg'],'jpg');
end

function plot_scatter(x1,x2,t)

k = (~isnan(x1)).*(~isnan(x2)) == 1;
mnX = min([x1(k);x2(k)]);
mxX = max([x1(k);x2(k)]);
hold on;
dscatter(x1,x2,'MSIZE',25);
line([mnX mxX],[mnX mxX],'LineStyle','-','color','k','linewidth',1);
hold off;
axis square;
set(gca,'xlim',[mnX mxX],'ylim',[mnX mxX],'fontsize',18);
xlabel('degradation model');
ylabel('polyA model');
title(sprintf('%s n=%d,r=%.2f',t,sum(k),corr(x1(k),x2(k))));
