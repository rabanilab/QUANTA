function h = fit_model_plot_param(P,Pids,PARAM_names,LB,UB)
% P = <logX0,dg,t0,t1>

n = max(size(P));
if (nargin < 3)
    PARAM_names = {'id' 'x0' 'dg' 't0' 't1'};
end
if (nargin < 5)
    LB = inf(1,4);
    UB = -inf(1,4);
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



h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(2,4,1);
plot_hist(LB(1),UB(1),P,Pids,1);
xlabel('logX0');

subplot(2,4,2);
plot_hist(LB(2),UB(2),P,Pids,2);
xlabel('dg rate (1/hr)');

subplot(2,4,6);
hlim = log(2)./[UB(2) LB(2)];
if (hlim(2)-hlim(1)<2)
    hlim = [0 10];
end
d = [];
for j = 1:n
    d{j} = log(2)./P{j}(:,2);
end
plot_hist(hlim(1),hlim(2),d,Pids,1);
xlabel('mRNA half-life (hr)');

subplot(2,4,3);
plot_hist(LB(3),UB(3),P,Pids,3);
xlabel('onset (hr)');

subplot(2,4,4);
plot_hist(LB(4),UB(4),P,Pids,4);
xlabel('offset (hr)');


function plot_hist(LB,UB,P,Pids,i)

if (UB-LB<1)
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
plot(x',y./sum(y),'-','linewidth',2);%,'marker','.','markersize',20);
axis tight;
%set(gca,'ylim',[0 maxy]);
ylabel('fraction');
box('off');
legend(L,'box','off','location','northOutside');
set(gca,'ytick',0:0.05:1,'fontsize',14);
