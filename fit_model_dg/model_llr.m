function [h,Sa,Pa,Ra,Ea] = model_llr(all_data_r, all_data_a, mtid, fpkm_pref)

llr_alpha = 0.05;
sig = 0.2;

% LLR test
nr = size(all_data_r,2);
na = size(all_data_a,2);
maxN = max(nr,na);

Sa = cell(1,nr+na);
Pa = cell(1,nr+na);
Ra = cell(1,nr+na);
Ea = cell(1,nr+na);
for j = 1:nr
    [mtid,Sa{j},Pa{j},Ra{j},Ea{j}] = model_select(mtid,fpkm_pref,all_data_r{j},llr_alpha,sig);
end
for j = 1:na
    [mtid,Sa{nr+j},Pa{nr+j},Ra{nr+j},Ea{nr+j}] = model_select(mtid,fpkm_pref,all_data_a{j},llr_alpha,sig);
end

% plot test results
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for j = 1:nr
    subplot(2,maxN,j);
    y = hist(Sa{j},0:2);
    b = bar(y./sum(y));
    axis tight;
    ylabel('fraction of genes');
    set(gca,'xticklabel',{'L0' 'L1' 'L2'},'ylim',[0 1],'fontsize',14);
    b.FaceColor = 'flat';
    for i = 1:3
        text(i-0.3,y(i)./sum(y)-0.05,sprintf('n=%d\n(%.1f%%)',y(i),100*y(i)./sum(y)),'fontsize',14);
        b.CData(i,:) = [0.8 0.8 0.8];
    end
    title(all_data_r{j});
end
for j = 1:na
    subplot(2,maxN,maxN+j);
    y = hist(Sa{nr+j},0:2);
    b = bar(y./sum(y));
    axis tight;
    ylabel('fraction of genes');
    set(gca,'xticklabel',{'L0' 'L1' 'L2'},'ylim',[0 1],'fontsize',14);
    b.FaceColor = 'flat';
    for i = 1:3
        text(i-0.3,y(i)./sum(y)-0.05,sprintf('n=%d\n(%.1f%%)',y(i),100*y(i)./sum(y)),'fontsize',14);
        b.CData(i,:) = [0.8 0.8 0.8];
    end
    title(all_data_a{j});
end


function [mtid,S,Pall,Rall,Eall] = model_select(mtid, fpkm_pref, data_id, llr_alpha, sig)

% expression data
[D,t,iM,~,~,tid] = load_data([fpkm_pref data_id '.txt']);
[max(size(t(iM))) min(t(iM)) max(t(iM))]
[~,i1,i2] = intersect(mtid,tid);
mtid = mtid(i1);
Y = D(i2(i2>0),iM);

% Null model [x0]
P0 = mean(Y,2);
X0 = repmat(P0,1,sum(iM));
R0 = zeros(size(P0));
E0 = sum((X0-Y).^2,2);
L0 = sum(log2(normpdf(Y,X0,sig)),2);

% Model #1 [x0,dg,t0,t1]
load(['model_param.1p.' data_id '.mat'],'mP','mR','mE');
P1 = mP(i1,:);
R1 = mR(i1,:);
E1 = mE(i1,:);
X1 = [];
for i = 1:size(P1,1)
    X1(i,:) = dg_eval_model(t(iM), P1(i,:));
end
L1 = sum(log2(normpdf(Y,X1,sig)),2); % ***
fprintf('rows with L0>L1: %d\n', sum(L1<L0))
[~,Q(:,1)] = lratiotest(L1,min(L0,L1),size(P1,2)-size(P0,2));
fprintf('rows reject L0 for L1: %d\n', sum(Q(:,1)<llr_alpha));

% Model #2 [x0,da,dg,t0,t1]
load(['model_param.2p.' data_id '.mat'],'mP','mR','mE');
P2 = mP(i1,:);
R2 = mR(i1,:);
E2 = mE(i1,:);
X2 = [];
for i = 1:size(P2,1)
    X2(i,:) = dg_eval_model_2p(t(iM), P2(i,:));
end
L2 = sum(log2(normpdf(Y,X2,sig)),2);
fprintf('rows with L1>L2: %d\n', sum(L2<L1));
[~,Q(:,2)] = lratiotest(L2,min(L1,L2),size(P2,2)-size(P1,2));
fprintf('rows reject L1 for L2: %d\n', sum(Q(:,2)<llr_alpha));

fprintf('rows with L0>L2: %d\n', sum(L2<L0));
[~,Q(:,3)] = lratiotest(L2,min(L0,L2),size(P2,2)-size(P0,2));
fprintf('rows reject L0 for L2: %d\n', sum(Q(:,3)<llr_alpha));

%figure;
%dscatter(-1*L1,-1*L2);
%set(gca,'xlim',[0 300],'ylim',[0 300]);
%axis square;
%xlabel('log likelihood M1');
%ylabel('log likelihood M2');

% select model per row
S = zeros(size(Q,1),1);
S(Q(:,1)<llr_alpha) = 1;
S((Q(:,1)<llr_alpha).*(Q(:,2)<llr_alpha)==1) = 2;
S((Q(:,1)>=llr_alpha).*(Q(:,3)<llr_alpha)==1) = 2;
y = hist(S,0:2);
num2cell([(0:2)' y' y'./sum(y)])

[u,~,t] = unique([Q<llr_alpha S],'rows');
c = accumarray(t,1);
num2cell([u c c./sum(c)])

Pall = zeros(size(P2));
Rall = zeros(size(R2));
Eall = zeros(size(E2));
k = (S==0);
Pall(k,:) = [P0(k) zeros(sum(k),4)];
Rall(k,:) = R0(k);
Eall(k,:) = E0(k);
k = (S==1);
Pall(k,:) = [P1(k,1) zeros(sum(k),1) P1(k,2:4)];
Rall(k,:) = R1(k);
Eall(k,:) = E1(k);
k = (S==2);
Pall(k,:) = P2(k,:);
Rall(k,:) = R2(k);
Eall(k,:) = E2(k);
