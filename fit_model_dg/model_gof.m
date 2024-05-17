function [h,p1r,p2r,p1a,p2a] = model_gof(all_data_r, all_data_a, mtid, fpkm_pref)

gof_alpha = 0.05;
sig = 0.2;

% GOF test
nmtid = size(mtid,1);
nr = size(all_data_r,2);
na = size(all_data_a,2);
[p1r,p2r] = model_compare(all_data_r,mtid,fpkm_pref,sig);
[p1a,p2a] = model_compare(all_data_a,mtid,fpkm_pref,sig);

% plot test results
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(1,2,1);
y = [sum(p1r<gof_alpha);sum(p2r<gof_alpha)]';
y = [y;sum(y,1)];
y = y./[repmat(nmtid,nr,1);nmtid*nr];
bar(y);
axis tight;
ylabel('fraction of genes rejecting the model');
if (max(y(:))<0.3)
    set(gca,'xticklabel',[all_data_r 'all'],'ylim',[0 0.3],'fontsize',14);
else
    set(gca,'xticklabel',[all_data_r 'all'],'ylim',[0 1],'fontsize',14);
end
for i = 1:nr
    text(i-0.3,y(i,1),sprintf('n=%d\n(%.1f%%)',round(nmtid*y(i,1)),100*y(i,1)),'fontsize',14);
    text(i-0.1,y(i,2),sprintf('n=%d\n(%.1f%%)',round(nmtid*y(i,2)),100*y(i,2)),'fontsize',14);
end
title('ribosomal depleted datasets');
legend({'dg model' 'pA model'},'location','northoutside','box','off');

subplot(1,2,2);
y = [sum(p1a<gof_alpha);sum(p2a<gof_alpha)]';
y = [y;sum(y,1)];
y = y./[repmat(nmtid,na,1);nmtid*na];
bar(y);
axis tight;
ylabel('fraction of genes rejecting the model');
if (max(y(:))<0.3)
    set(gca,'xticklabel',[all_data_a 'all'],'ylim',[0 0.3],'fontsize',14);
else
    set(gca,'xticklabel',[all_data_a 'all'],'ylim',[0 1],'fontsize',14);
end
for i = 1:na
    text(i-0.3,y(i,1),sprintf('n=%d\n(%.1f%%)',round(nmtid*y(i,1)),100*y(i,1)),'fontsize',14);
    text(i-0.1,y(i,2),sprintf('n=%d\n(%.1f%%)',round(nmtid*y(i,2)),100*y(i,2)),'fontsize',14);
end
title('polyA selected datasets');
legend({'dg model' 'pA model'},'location','northoutside','box','off');



function [p1,p2] = model_compare(all_data, mtid, fpkm_pref, sig)

n = size(all_data,2);

p1 = [];
s = [];
for j = 1:n
    [D,t,iM,~,~,tid] = load_data([fpkm_pref all_data{j} '.txt']);
    load(['model_param.1p.' all_data{j} '.mat'],'mE');
    i = ismember(tid,mtid);
    W = D(i,iM);
    s(j) = sig*std(W(:),1).^2;
    k = size(t,2);
    p1(:,j) = chi2cdf(mE./s(j),k,'upper');
end

p2 = [];
s = [];
for j = 1:n
    [D,t,iM,~,~,tid] = load_data([fpkm_pref all_data{j} '.txt']);
    load(['model_param.2p.' all_data{j} '.mat'],'mE');
    i = ismember(tid,mtid);
    W = D(i,iM);
    s(j) = sig*std(W(:),1).^2;
    k = size(t,2);
    p2(:,j) = chi2cdf(mE./s(j),k,'upper');
end
    
