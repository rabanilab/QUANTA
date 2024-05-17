function h = plot_compare_class(n,orgid,CLASS,orgid2)

if (nargin < 4)
    orgid2 = 'zfish';
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
c = [1 0 0; 0 0 1;0.8 0.8 0.8];

subplot(1,2,1);
b = bar(n./sum(n(:)),'stacked');
for i = 1:3
    b(i).FaceColor = 'flat';
    p = cumsum(n(i,:));
    for j = 1:3
        text(i-0.3,p(j)./sum(n(:)),sprintf('n=%d (%.1f%%)',n(i,j),100*n(i,j)./sum(n(i,:))),'fontsize',14);
        b(j).CData(i,:) = c(j,:);
    end
end
axis tight;
set(gca,'xticklabel',CLASS,'ylim',[0 1],'fontsize',18);
ylabel('fraction of genes');
title([orgid2 ' classification']);
legend(strcat(orgid,CLASS),'box','off');

subplot(1,2,2);
b = bar(n'./sum(n(:)),'stacked');
for i = 1:3
    b(i).FaceColor = 'flat';
    p = cumsum(n(:,i));
    for j = 1:3
        text(i-0.3,p(j)./sum(n(:)),sprintf('n=%d (%.1f%%)',n(j,i),100*n(j,i)./sum(n(:,i))),'fontsize',14);
        b(j).CData(i,:) = c(j,:);
    end
end
axis tight;
set(gca,'xticklabel',CLASS,'ylim',[0 1],'fontsize',18);
ylabel('fraction of genes');
title(sprintf('%s classification',orgid));
legend(strcat(orgid2,CLASS),'box','off');
