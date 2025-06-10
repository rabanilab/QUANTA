function plot_known_kmers(Sid,S,Kids,D,Dids,output_id,Dlab,Dlim,Dres,plotType,plotEPS)

if (nargin < 7)
    Dlab = 'avg. polyA length';
end
if (nargin < 8)
    maxD = nanmax(D(:));
    minD = nanmin(D(:));
    Dlim = [minD maxD];
end
if (nargin < 9)
    Dres = 2;
end
if (nargin < 10)
    plotType = 1;
end
if (nargin < 11)
    plotEPS = 0;
end

len_fold = 10;

% find kmer occurance
Slen = cellfun(@length,S);
n = max(size(Kids));
K = [];
for i = 1:n
    K(:,i) = cellfun(@isempty,regexp(S,Kids{i})) == 0;
end
oids = strrep(strrep(Kids,'(',''),')','');
oids = strrep(strrep(oids,'{',''),'}','');
oids = strrep(strrep(oids,'[',''),']','');
oids = strrep(oids,'|','.');
oids = strrep(oids,'*','.');

I = (Slen > (1/len_fold)*median(Slen)).*(Slen < len_fold*median(Slen)) == 1;
K = K(I,:);
D = D(I,:);

Kx = [];
if (plotType == 2)
    for i = 1:n
        Kx(:,i) = K(:,i).*(sum(K(:,setdiff(1:3,i)),2)==0) == 1;
    end
    [Kids' num2cell([sum(K);sum(K==0);sum(Kx)]')]
else
    [Kids' num2cell([sum(K);sum(K==0)]')]
end

% make plots
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:size(K,2)
    if (plotType == 2)
        plot_hist(Kids{i},Kx(:,i),K(:,i)==0,D,Dids,Dlim,Dlab,Dres);
    elseif (plotType == 3)
        plot_box(Kids{i},K(:,i),D,Dids,Dlim,Dlab,1,0);
    elseif (plotType == 4)
        plot_box(Kids{i},K(:,i),D,Dids,Dlim,Dlab,0,1);
    elseif (plotType == 5)
        plot_violin(Kids{i},K(:,i),D,Dids,Dlim,Dlab);
    else
        plot_hist(Kids{i},K(:,i),K(:,i)==0,D,Dids,Dlim,Dlab,Dres);
    end

    if (plotEPS > 0)
        saveas(h, [output_id '.' oids{i} '.svg'],'svg');
    else
        saveas(h, [output_id '.' oids{i} '.jpg'],'jpg');
    end
end
close all;


function plot_box(Kid,K,D,Dids,Dlim,Dlab,plotLine,plotBar)

clf;
n = size(D,2);
Dids = regexprep(Dids,'_',' ');

if (plotBar)
    hold on;
    boxplot(D(K==0,:),'boxstyle','outline','Symbol','.','widths',0.3,'Colors','k','positions',(1:n) + 0.35);
    boxplot(D(K==1,:),'boxstyle','outline','Symbol','.','widths',0.3,'Colors','k','positions',1:n,'labels',Dids);
    q = findobj(gca,'Tag','Box');
    for i = 1:size(D,2)
        patch(get(q(i),'XData'),get(q(i),'YData'),[1 0 0],'FaceAlpha',.5);
    end
    for i = size(D,2)+1:length(q)
        patch(get(q(i),'XData'),get(q(i),'YData'),[0.2 0.2 0.2],'FaceAlpha',.5);
    end
    hold off;
    q = findobj(gca,'type','line');
    set(q,'linew',1);
    q = findobj('LineStyle','--');
    set(q,'LineStyle','-');
end
if (plotLine)
    hold on;
    x = (1:n) + 0.35/2;
    m = nanmedian(D(K==0,:));
    plot(x,m,'-k','linewidth',2);
    m = nanmedian(D(K==1,:));
    plot(x,m,'-r','linewidth',2);
    hold off;
    set(gca,'ylim',Dlim);
end
ylabel(Dlab);
set(gca,'fontsize',15,'xlim',[0.75 n+0.75]);%,'ylim',Dlim
box off;
for i = 1:n
    if ((sum(K==0)>0)*(sum(K==1)>0)==1)
        e = nanmean(D(K==1,i)) - nanmean(D(K==0,i));
        if (e<0)
            [~,q1] = kstest2(D(K==0,i),D(K==1,i),'tail','smaller');
            y = Dlim(1) + 0.2*(Dlim(2)-Dlim(1));
            text(i,y,{sprintf('%.0e',q1) sprintf('%.1f',e)},'FontSize',15);
        elseif (e>0)
            [~,q2] = kstest2(D(K==0,i),D(K==1,i),'tail','larger');
            y = Dlim(1) + 0.8*(Dlim(2)-Dlim(1));
            text(i,y,{sprintf('%.0e',q2) sprintf('%.1f',e)},'FontSize',15);
        end
    end
end
title(sprintf('%s (n=%d, n0=%d)',Kid,sum(K==1),sum(K==0)))


function plot_violin(Kid,K,D,Dids,Dlim,Dlab)

clf;
n = size(D,2);
Dids = regexprep(Dids,'_',' ');

hold on;
w = [1:n;n+1:2*n];
x = [D(K==0,:) zeros(size(D(K==0,:)))];
x = x(:,w(:));
y = repmat(1:2*n,sum(K==0),1);
violinplot(x(:),y(:),'ShowData',false,'ViolinAlpha',1,'ViolinColor',[0.9 0.1 0.1],'Bandwidth',0.1);
x = [zeros(size(D(K==1,:))) D(K==1,:)];
x = x(:,w(:));
y = repmat(1:2*n,sum(K==1),1);
violinplot(x(:),y(:),'ShowData',false,'ViolinAlpha',1,'ViolinColor',[0.9 0.9 0.9],'Bandwidth',0.1);
hold off;
axis tight;

ylabel(Dlab);
set(gca,'fontsize',15);
set(gca,'xtick',(1:2:(2*n-1))+0.5,'xticklabel',Dids);
box off;
for i = 1:n
    if ((sum(K==0)>0)*(sum(K==1)>0)==1)
        e = nanmean(D(K==1,i)) - nanmean(D(K==0,i));
        if (e<0)
            [~,q1] = kstest2(D(K==0,i),D(K==1,i),'tail','smaller');
            y = Dlim(1) + 0.2*(Dlim(2)-Dlim(1));
            text(2*i-1,y,{sprintf('%.0e',q1) sprintf('%.1f',e)},'FontSize',15);
        elseif (e>0)
            [~,q2] = kstest2(D(K==0,i),D(K==1,i),'tail','larger');
            y = Dlim(1) + 0.8*(Dlim(2)-Dlim(1));
            text(2*i-1,y,{sprintf('%.0e',q2) sprintf('%.1f',e)},'FontSize',15);
        end
    end
end
title(sprintf('%s (n=%d, n0=%d)',Kid,sum(K==1),sum(K==0)))


function plot_hist(Kid,K1,K0,D,Did,Dlim,Dlab,Dres)

clf;
Did = regexprep(Did,'_',' ');
x = Dlim(1):Dres:Dlim(2);
n = min(size(D,2),3*4);
for i = 1:n
    subplot(3,4,i);
    %subplot(2,2,i);
    i1 = K1.*(isnan(D(:,i))==0)==1;
    i0 = K0.*(isnan(D(:,i))==0)==1;
    m1 = 0;
    m0 = 0;
    clear y;
    hold on;
    if (sum(i0)>0)
        y = hist(D(i0,i),x);
        plot(x,y./sum(y,2),':k','marker','.','markersize',20,'linewidth',2);
        m0 = nanmean(D(i0,i));
    end
    if (sum(i1)>0)
        y = hist(D(i1,i),x);
        plot(x,y./sum(y,2),'-r','marker','.','markersize',20,'linewidth',2);
        m1 = nanmean(D(i1,i));
    end
    hold off;
    axis tight;
    xlabel(Dlab);
    legend({sprintf('other (%.1f, %d)',m0,sum(i0)), ...
        sprintf('+seq (%.1f, %d)',m1,sum(i1))}, ...
        'box','off');
    if ((sum(i1)>0)*(sum(i0)>0)==1)
        e = m1 - m0;
        if (e<0)
            [~,q1] = kstest2(D(i0,i),D(i1,i),'tail','smaller');
            title(sprintf('%s (p=%.0e, e=%.1f, smaller)',Did{i},q1,e));
        elseif (e>0)
            [~,q2] = kstest2(D(i0,i),D(i1,i),'tail','larger');
            title(sprintf('%s (p=%.0e, e=%.1f, larger)',Did{i},q2,e));
        else
            title(sprintf('%s',Did{i}));
        end
    else
        title(sprintf('%s',Did{i}));
    end
    set(gca,'fontsize',15);
end
