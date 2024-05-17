function h = fit_model_2p_plot_param(P,Pids,PARAM_names,LB,UB)
% P = <logX0,dg1,dg2,t0>

n = max(size(P));
if (nargin < 3)
    PARAM_names = {'id' 'x0' 'dg' 'da' 't0' 't1'};
end
if (nargin < 5)
    LB = inf(1,5);
    UB = -inf(1,5);
    for j = 1:n
        mx = max(P{j});
        mn = min(P{j});
        LB = min(LB,mn);
        UB = max(UB,mx);
    end
end


fprintf('MEAN PARAM:\n');
W = [];
for j = 1:n
    W = [W;[Pids(j) num2cell(nanmean(P{j}))]];
end
[PARAM_names;W]

fprintf('3-th prctile PARAM:\n');
W = [];
for j = 1:n
    W = [W;[Pids(j) num2cell(prctile(P{j},3))]];
end
[PARAM_names;W]

fprintf('97-th prctile PARAM:\n');
W = [];
for j = 1:n
    W = [W;[Pids(j) num2cell(prctile(P{j},97))]];
end
[PARAM_names;W]

fprintf('MEDIAN PARAM:\n');
W = [];
for j = 1:n
    W = [W;[Pids(j) num2cell(nanmedian(P{j}))]];
end
[PARAM_names;W]

H = [];
for j = 1:n
    fprintf('MEDIAN half-life %s:\n',Pids{j});
    hl = 60*log(2)./P{j}(:,3);
    H = [H;hl];
    [nanmedian(hl) num2cell(nanmedian(H(:)))]
end


h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
subplot(2,5,1);
plot_hist(LB(1),UB(1),P,Pids,1);
xlabel('logX0');

subplot(2,5,2);
plot_hist(LB(2),UB(2),P,Pids,2);
xlabel('da rate (1/hr)');

subplot(2,5,7);
d = [];
mnd = 1;
for j = 1:n
    prctile(abs(P{j}(:,2)),20)
    d{j} = log(2)./P{j}(:,2);
    mnd = min(mnd,prctile(abs(P{j}(:,2)),20));
end
if (mnd == 0)
    mnd = 0.01;
end
hlim = log(2)./[-1*mnd mnd];
if (hlim(2)-hlim(1)<2)
    hlim = [-10 10];
end
plot_hist(hlim(1),hlim(2),d,Pids,1);
xlabel('A+ half-life (hr)');

subplot(2,5,3);
plot_hist(LB(3),UB(3),P,Pids,3);
xlabel('dg rate (1/hr)');

subplot(2,5,8);
hlim = log(2)./[UB(3) LB(3)];
if (hlim(2)-hlim(1)<2)
    hlim = [0 10];
end
d = [];
for j = 1:n
    d{j} = log(2)./P{j}(:,3);
end
plot_hist(hlim(1),hlim(2),d,Pids,1);
xlabel('mRNA half-life (hr)');

subplot(2,5,4);
plot_hist(LB(4),UB(4),P,Pids,4);
xlabel('onset (hr)');

subplot(2,5,5);
plot_hist(LB(5),UB(5),P,Pids,5);
xlabel('offset (hr)');


function plot_hist(LB,UB,P,Pids,i)

if (UB-LB<0.5)
    UB = UB+1;
end
hd = round(abs(LB-UB)/35,2);
x = LB:hd:UB;

n = size(P,2);
L = [];
for j = 1:n
    d = P{j}(:,i);
    k = (~isnan(d)).*(~isinf(d)).*(d~=0) == 1;
    y(:,j) = hist(d(k),x)';
    L{j} = sprintf('%s (n=%d,m=%.2f)', regexprep(Pids{j},'_',' '), sum(k),mean(d(k)));
end
plot(x',y./sum(y),'-','linewidth',2);
axis tight;
ylabel('fraction');
box('off');
legend(L,'box','off','location','northOutside');
if (max(max(y./sum(y)))>0.2)
    set(gca,'ytick',0:0.1:1,'fontsize',14);
else
    set(gca,'ytick',0:0.05:1,'fontsize',14);
end
