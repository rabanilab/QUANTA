function h = plot_ATR(xT,Xa,Xr,Rx)

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(1,2,1);
hold on;
plot(xT,mean(Xa),'-r','marker','.','LineWidth',2,'markersize',20);
plot(xT,mean(Xr),'-b','marker','.','LineWidth',2,'markersize',20);
hold off;
xlabel('time (hpf)');
ylabel('model predicted expression (mean)');
legend({'polyA' 'total'},'location','bestOutside','box','off');
set(gca,'fontsize',18);

subplot(1,2,2);
hold on;
plot(xT,mean(Rx),'-k','marker','.','LineWidth',2,'markersize',20);
line([min(xT) max(xT)],[0 0],'LineStyle','-','color','r','linewidth',1.2);
hold off;
xlabel('time (hpf)');
ylabel('expression ratio (mean)');
set(gca,'ylim',[-2 2],'fontsize',18);
