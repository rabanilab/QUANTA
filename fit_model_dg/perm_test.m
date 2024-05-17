function [pv,v0,v1,h] = perm_test(C,n)

if (nargin < 2)
    n = 10000;
end
m = size(C,1);

% background distribution by random permutations
v0 = [];
for i = 1:n
    j = randperm(m)';
    for w = 1:3
        v0(i,w) = sum((C(:,1)==w).*(C(j,2)==w));
    end
    v0(i,4) = sum(C(:,1)==C(j,2));
end

% actual counts
v1 = [];
for w = 1:3
    v1(w) = sum((C(:,1)==w).*(C(:,2)==w));
end
v1(4) = sum(C(:,1)==C(:,2));

% pvalues (using normal distribution)
for w = 1:3
    pv(w) = normcdf(v1(w),mean(v0(:,w)),std(v0(:,w)),'upper');
end
pv(4) = normcdf(v1(4),mean(v0(:,4)),std(v0(:,4)),'upper');

% plot the results
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for w = 1:3
    subplot(2,2,w);
    mn = min(v0(:,w),v1(w));
    mx = max(v0(:,w),v1(w));
    x = max(mn-1000,0):5:min(mx+1000,m);
    hold on;
    y = hist(v0(:,w),x);
    plot(x,y./sum(y),'-k','linewidth',1.2);
    line([v1(w) v1(w)],[0 max(y)/sum(y)],'LineStyle','-','color','r','linewidth',2);
    hold off;
    xlabel('matched labels');
    ylabel('fraction');
    legend({'random tests' 'actual'},'Location','bestoutside','box','off');
    set(gca,'fontsize',18);
    title(sprintf('(%d) p<%.2e',w,pv(w)));
end
subplot(2,2,4);
mn = min(v0(:,4),v1(4));
mx = max(v0(:,4),v1(4));
x = max(mn-100,0):5:min(mx+100,m);
hold on;
y = hist(v0(:,4),x);
plot(x,y./sum(y),'-k','linewidth',1.2);
line([v1(4) v1(4)],[0 max(y)/sum(y)],'LineStyle','-','color','r','linewidth',2);
hold off;
xlabel('matched labels');
ylabel('fraction');
legend({'random tests' 'actual'},'Location','bestoutside','box','off');
set(gca,'fontsize',18);
title(sprintf('p<%.2e',pv(4)));
