function h = compare_pa_rz(mtid,all_data_r,Pr,all_data_a,Pa)
% compare ribo to polya

mN = size(mtid,1);

% average parameters across all datasets
mPr = zeros(mN,size(Pr{1},2));
Cr = zeros(mN,size(Pr{1},2));
for j = 1:max(size(all_data_r))
    for i = 1:size(Pr{j},2)
        mPr(:,i) = mPr(:,i) + Pr{j}(:,i);
        k = (~isnan(Pr{j}(:,i))).*(~isinf(Pr{j}(:,i))).*(Pr{j}(:,i) ~= 0);
        Cr(:,i) = Cr(:,i) + k;
    end
end
mPr = mPr./Cr;

mPa = zeros(mN,size(Pa{1},2));
Ca = zeros(mN,size(Pa{1},2));
for j = 1:max(size(all_data_a))
    for i = 1:size(Pa{j},2)
        mPa(:,i) = mPa(:,i) + Pa{j}(:,i);
        k = (~isnan(Pa{j}(:,i))).*(~isinf(Pa{j}(:,i))).*(Pa{j}(:,i) ~= 0);
        Ca(:,i) = Ca(:,i) + k;
    end
end
mPa = mPa./Ca;

% plot comparisons
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

x1 = mPr(:,1);
x2 = mPa(:,1);
plot_corr(x1,'x0 ribo',x2,'x0 polya',1);
x1 = log2(mPr(:,3));
x2 = log2(mPa(:,3));
plot_corr(x1,'dg ribo; log2',x2,'dg polya; log2',2);
x1 = mPr(:,4);
x2 = mPa(:,4);
plot_corr(x1,'t0 ribo',x2,'t0 polya',3);

x1 = log2(mPr(:,3));
x2 = mPa(:,1) - mPr(:,1);
plot_corr(x1,'dg (ribo, log2)',x2,'x0 ratio',5);
axis tight;
x2 = mPa(:,2);
plot_corr(x1,'dg (ribo, log2)',x2,'da polya',6);
axis tight;
x1 = log2(mPa(:,3));
x2 = mPa(:,1) - mPr(:,1);
plot_corr(x1,'dg (polyA, log2)',x2,'x0 ratio',7);
axis tight;
x2 = mPa(:,2);
plot_corr(x1,'dg (polyA, log2)',x2,'da polya',8);
axis tight;


function plot_corr(x1,t1,x2,t2,j,maxN)

if (nargin < 6)
    maxN = 4;
end

subplot(2,maxN,j);
k = (~isnan(x1)).*(~isnan(x2)).*(~isinf(x1)).*(~isinf(x2))==1;
if (sum(k) > 10)
    rn1 = 0.01*randn(size(x1(k)));
    rn2 = 0.01*randn(size(x2(k)));
    dscatter(x1(k)+rn1,x2(k)+rn2,'MSIZE',20);
    minX = min([x1(k);x2(k)]);
    maxX = max([x1(k);x2(k)]);
    set(gca,'xlim',[minX maxX],'ylim',[minX maxX]);
    r = unique(round(minX:maxX));
    if (max(size(r))>5)
        dr = round(max(size(r))/5);
        if (dr<1)
            dr = 1;
        end
        r = r(1:dr:end);
    end
    set(gca,'xtick',r,'ytick',r,'fontsize',14);
    title(sprintf('r=%.2f, n=%d',corr(x1(k),x2(k)),sum(k)));
end
xlabel(regexprep(t1,'_',' '));
ylabel(regexprep(t2,'_',' '));
axis square;

