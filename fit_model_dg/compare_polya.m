function compare_polya(Mid,M,t,mtid,Ma,La,out_prefix,ylim)

nm = size(M,2);
nma = size(Ma,2);
[~,j1,j2] = intersect(Mid,mtid);
sum(strcmp(Mid(j1),mtid(j2(j2>0))))

% correlation
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(2,4,1);
c = corr(Ma(j2(j2>0),:),log2(M(j1,:)),'rows','pairwise');
imagesc(c,[-1 1]);
colormap(gene_colormap(1));
axis square;
set(gca,'xtick',1:nm,'xticklabel',t,'ytick',1:nma,'yticklabel',La);
xlabel('pA lengths');
ylabel('model');
colorbar;
for i = 1:nm
    subplot(2,4,i+1);
    dM = log2(M) - log2(M(:,i));
    c = corr(Ma(j2(j2>0),:),dM(j1,:),'rows','pairwise');
    imagesc(c,[-1 1]);
    colormap(gene_colormap(1));
    axis square;
    set(gca,'xtick',1:nm,'xticklabel',t,'ytick',1:nma,'yticklabel',La);
    xlabel(sprintf('pA diff %d',t(i)));
    ylabel('model');
    colorbar;
end
saveas(h,[out_prefix '.corr.jpg'],'jpg');

colormap default;
for k = 1:nm
    clf;
    for i = 1:nma
        subplot(3,4,i);
        i1 = log2(M(j1,k));
        i2 = Ma(j2(j2>0),i);
        j = (~isnan(i1)).*(~isnan(i2))==1;
        dscatter(i1(j),i2(j),'MSIZE',20);
        axis square;
        xlabel(sprintf('pA %d hr; log2', t(k)));
        ylabel(La{i});
        set(gca,'xlim',[0 8],'ylim',ylim);
        set(gca,'fontsize',14);
        title(sprintf('(r=%.2f,n=%d)',corr(i1(j),i2(j)),sum(j)));
    end
    saveas(h,[out_prefix '.' num2str(t(k)) 'h.jpg'],'jpg');
end
for k = 2:nm
    clf;
    for i = 1:nma
        subplot(3,4,i);
        i1 = log2(M(j1,1))-log2(M(j1,k));
        i2 = Ma(j2(j2>0),i);
        j = (~isnan(i1)).*(~isnan(i2))==1;
        dscatter(i1(j),i2(j),'MSIZE',20);
        axis square;
        xlabel(sprintf('pA %d hr - %d hr; log2', t(k),t(1)));
        ylabel(La{i});
        set(gca,'xlim',[-5 5],'ylim',ylim);
        set(gca,'fontsize',14);
        title(sprintf('(r=%.2f,n=%d)',corr(i1(j),i2(j)),sum(j)));
    end
    saveas(h,[out_prefix '.diff.' num2str(t(k)) 'h.jpg'],'jpg');
end
close all;

