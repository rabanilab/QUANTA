function plot_known_kmers(Sid,S,Kids,D,Dids,output_id,Dlab,Dlim,Dres,plotType)

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

% find kmer occurance
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

% plot by kmer
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

if (plotType == 2)
    Kx = [];
    for i = 1:n
        Kx(:,i) = K(:,i).*(sum(K(:,setdiff(1:3,i)),2)==0) == 1;
    end
    [Kids' num2cell([sum(K);sum(K==0);sum(Kx)]')]
    for i = 1:size(K,2)
        Kids(i)
        plot_distributions(Kx(:,i),K(:,i)==0,D,Dids,Dlim,Dlab,Dres);
        saveas(h, [output_id '.' oids{i} '.jpg'],'jpg');
    end
elseif (plotType == 3)
    for i = 1:size(K,2)
        Kids(i)
        plot_boxplot(K(:,i),D,Dids,Dlim,Dlab,1);
        saveas(h, [output_id '.' oids{i} '.jpg'],'jpg');
    end
elseif (plotType == 4)
    for i = 1:size(K,2)
        Kids(i)
        plot_boxplot(K(:,i),D,Dids,Dlim,Dlab,0);
        saveas(h, [output_id '.' oids{i} '.jpg'],'jpg');
    end
else
    for i = 1:size(K,2)
        Kids(i)
        plot_distributions(K(:,i),K(:,i)==0,D,Dids,Dlim,Dlab,Dres);
        saveas(h, [output_id '.' oids{i} '.jpg'],'jpg');
    end
end

close all;


function plot_boxplot(K,D,Dids,Dlim,Dlab,plotLine)

clf;
Dids = regexprep(Dids,'_',' ');
hold on;
n = size(D,2);
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
if (plotLine)
    hold on;
    x = (1:n) + 0.35/2;
    m = nanmedian(D(K==0,:));
    plot(x,m,'-k','linewidth',2);
    m = nanmedian(D(K==1,:));
    plot(x,m,'-r','linewidth',2);
    hold off;
end
ylabel(Dlab);
set(gca,'fontsize',15);%'ylim',Dlim);
box off;


function plot_distributions(K1,K0,D,Did,Dlim,Dlab,Dres)

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
        m0 = nanmedian(D(i0,i));
    end
    if (sum(i1)>0)
        y = hist(D(i1,i),x);
        plot(x,y./sum(y,2),'-r','marker','.','markersize',20,'linewidth',2);
        m1 = nanmedian(D(i1,i));
    end
    hold off;
    axis tight;
    xlabel(Dlab);
    legend({sprintf('other (%.1f, %d)',m0,sum(i0)), ...
        sprintf('+seq (%.1f, %d)',m1,sum(i1))}, ...
        'box','off');
    if ((sum(i1)>0)*(sum(i0)>0)==1)
        [~,q1] = kstest2(D(i0,i),D(i1,i),'tail','smaller');
        [~,q2] = kstest2(D(i0,i),D(i1,i),'tail','larger');
        if (q1<q2)
            title(sprintf('%s (p=%.1e, smaller)',Did{i},q1));
        else
            title(sprintf('%s (p=%.1e, larger)',Did{i},q2));
        end
    else
        title(sprintf('%s',Did{i}));
    end
    set(gca,'fontsize',15);
end
