function compare_to_zebrafish(orgid,mzt,drange,all_names)

if (nargin < 3)
    drange = [-7 2];
end
if (nargin < 4)
    all_names = {'dazl' 'sod2' 'h1m' 'btg4' 'cpeb1' 'zar1' 'buc' 'wee2'};
end

zfish_mzt = 3.5;

% compare class
system('make compare_class_to_zfish');
X = importdata('orthologies.class.txt');

C = X.data;
write_text_file('zf_comp.class.txt',[X.textdata num2cell(C(:,2)+0.1*C(:,1))]);

[~,~,~,h] = perm_test(C);
saveas(h,'results/orthologies.class.perm.jpg','jpg');
saveas(h,'results/orthologies.class.perm.eps','epsc');

fprintf('%s:\n',orgid);
y = hist(C(:,1),1:3);
num2cell([y' y'./sum(y)])
fprintf('fish:\n');
y = hist(C(:,2),1:3);
num2cell([y' y'./sum(y)])
n = [];
for i = 1:3
    [u,~,t] = unique(C(C(:,2)==i,:),'rows');
    n(i,:) = accumarray(t,1);
    num2cell([u n(i,:)' n(i,:)'./sum(n(i,:))])
end

CLASS = {'M' 'Z' 'MZ'};
h = plot_compare_class(n,orgid,CLASS);
saveas(h,'results/orthologies.class.jpg','jpg');
saveas(h,'results/orthologies.class.eps','epsc');


% compare all rates
system('make compare_param_to_zfish');
Y = importdata('orthologies.x0.txt');
nx = Y.textdata(:,1);
%nx = Y.textdata(sum(Y.data>2,2)==2,1);

h = plot_compare_rates(nx,orgid,zfish_mzt/mzt,drange);
saveas(h, 'results/orthologies.jpg','jpg');

% compare degradation rates
X = importdata('orthologies.dg.txt');
i = ismember(X.textdata(:,1),nx);
x = log2(X.data(i,1));
y = log2(X.data(i,2));
id = X.textdata(i,:);
j = (~isnan(x)).*(~isnan(y)) == 1;
x = x(j);
y = y(j);

f = mean(y-x);
p = prctile(abs(x-(y-f)),25);
d = (x-(y-f) > p) + 2*((y-f)-x > p); % 0 = same; 1 = frog > fish; 2 = fish > frog
write_text_file('zf_comp.dg.txt',[id(j,:) num2cell(d)]);

A = importdata('../1_zebrafish_Nov2023/gene_classification.txt');
[~,j2] = ismember(id(j,2),A.textdata(:,1));
zid = A.textdata(j2(j2>0),2);

h = plot_compare_dg(x,y,zid,all_names,orgid,drange);
saveas(h, 'results/orthologies.dg.svg','svg');

close all;
