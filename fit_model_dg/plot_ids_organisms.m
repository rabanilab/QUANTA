function plot_ids_organisms(ids,data_a,data_r,loadpref,mzt,outpref,max_scaleT,scale_time,run_mfit,gene_ids,feps)

if (nargin < 7)
    max_scaleT = 10;
end
if (nargin < 8)
    scale_time = 0;
end
if (nargin < 9)
    run_mfit = 1;
end
if (nargin < 10)
    gene_ids = 0;
end
if (nargin < 11)
    feps = 0;
end

N = max(size(data_a));
y_scale = [-3 12];

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:max(size(ids))
    clf;
    C = colormap('lines');

    % polyA datasets
    subplot(1,2,1);
    mxT = 1;
    Q = [];
    S = [];
    T = [];
    hold on;
    for j = 1:N
        [gname1,Q1,S1,mxTj] = plot_a(ids{i,j},data_a{j},loadpref{j},scale_time,max_scaleT,run_mfit,gene_ids,mzt{j},C(j,:));
        mxT = max(mxT,mxTj);
        Q = [Q Q1];
        S = [S S1];
        if (isempty(T))
            T = gname1;
        else
            T = [T ',' gname1];
        end
    end
    hold off;
    if (scale_time)
        xlabel('scaled time (a.u.)');
    else
        xlabel('time (hr)');
    end
    ylabel('logFPKM');
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([T ' (polyA)']);
    set(gca,'ylim',y_scale,'xlim',[0 mxT+0.5],'fontsize',15);
    axis square;

    % ribo depletion datasets
    subplot(1,2,2);
    mxT = 1;
    Q = [];
    S = [];
    T = [];
    hold on;
    for j = 1:N
        [gname1,Q1,S1,mxTj] = plot_r(ids{i,j},data_r{j},loadpref{j},scale_time,max_scaleT,run_mfit,gene_ids,mzt{j},C(j,:));
        mxT = max(mxT,mxTj);
        Q = [Q Q1];
        S = [S S1];
        if (isempty(T))
            T = gname1;
        else
            T = [T ',' gname1];
        end
    end
    hold off;
    if (scale_time)
        xlabel('scaled time (a.u.)');
    else
        xlabel('time (hr)');
    end
    ylabel('logFPKM');
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([T ' (ribo)']);
    set(gca,'ylim',y_scale,'xlim',[0 mxT+0.5],'fontsize',15);
    axis square;

    % save file
    if (~isempty(gname1))
        saveas(h, [outpref gname1 '.jpg'],'jpg');
        if (feps > 0)
            saveas(h, [outpref gname1 '.svg'],'svg');
        end
    else
        saveas(h, [outpref ids{i} '.EMPTY.jpg'],'jpg');
    end
end

close all;


function [gname,Q,S,mxT] = plot_a(id,data_a,loadpref,scale_time,max_scaleT,run_mfit,gene_ids,mzt,plotcolor)

name_a = regexprep(data_a,'_',' ');

maxT = 1;
if (scale_time)
    for j = 1:max(size(data_a))
        [~,t] = load_data([loadpref data_a{j} '.txt']);
        maxT = max(maxT,max(t));
    end
end

S = {};
Q = [];
gname = [];
mxT = 1;
for j = 1:max(size(data_a))
    [D,t,iM,iP,gid,tid,gc] = load_data([loadpref data_a{j} '.txt']);
    mxT = max(mxT,max((max_scaleT/maxT).*t));
    
    if (gene_ids > 0)
        k = find(strcmp(gid,id));
    else
        k = find(strcmp(tid,id));
    end
    if ((~isempty(k))&&(max(size(k))>1))
        fprintf('** gene %s found %d lines **\n', id, max(size(k)));
        k = k(1);
    end

    if ((~isempty(k))&&(max(size(k))==1))
        gname = gid{k};
        x = t(iM);
        tX = min(x):0.1:max(x);
        y = D(k,iM);
        plot((max_scaleT/maxT).*x,y,'.','markersize',35,'color',plotcolor);
        if (run_mfit)
            if (size(x,2)>=5)
                P = fit_model_2p(x,y,mzt);
            elseif (size(x,2)>=4)
                P = fit_model(x,y,mzt(2:3));
                P = [P(1) 0 P(2:4)];
            else
                P = [mean(y) 0 0 0 0];
            end
            X = dg_eval_model_2p(tX,P);
            Q(j) = plot((max_scaleT/maxT).*tX,X,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s: %.1f min, %.1f min, %.1f hr, %.1f]',...
                name_a{j},gc{k},60*log(2)/P(2),60*log(2)/P(3),P([4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot((max_scaleT/maxT).*u,s,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s]',name_a{j},gc{k});
        end
        %z = D(k,iP);
        %plot(x./mxS,z,'s','markersize',10,'color',plotcolor);
        %[u,~,w] = unique(x);
        %s = accumarray(w,z)./accumarray(w,1);
        %plot(u./mxS,s,'--','linewidth',2,'color',plotcolor);
    end
end


function [gname,Q,S,mxT] = plot_r(id,data_r,loadpref,scale_time,max_scaleT,run_mfit,gene_ids,mzt,plotcolor)

name_r = regexprep(data_r,'_',' ');

maxT = 1;
if (scale_time)
    for j = 1:max(size(data_r))
        [~,t] = load_data([loadpref data_r{j} '.txt']);
        maxT = max(maxT,max(t));
    end
end

S = {};
Q = [];
gname = [];
mxT = 1;
for j = 1:max(size(data_r))
    [D,t,iM,iP,gid,tid,gc] = load_data([loadpref data_r{j} '.txt']);
    mxT = max(mxT,max((max_scaleT/maxT).*t));

    if (gene_ids > 0)
        k = find(strcmp(gid,id));
    else
        k = find(strcmp(tid,id));
    end
    if ((~isempty(k))&&(max(size(k))>1))
        fprintf('** gene %s found %d lines **\n', id, max(size(k)));
        k = k(1);
    end

    if ((~isempty(k))&&(max(size(k))==1))
        gname = gid{k};
        x = t(iM);
        tX = min(x):0.1:max(x);
        y = D(k,iM);
        plot((max_scaleT/maxT).*x,y,'.','markersize',35,'color',plotcolor);
        if (run_mfit)
            if (size(x,2)>=4)
                P = fit_model(x,y,mzt(2:3));
            else
                P = [mean(y) 0 0 0];
            end
            X = dg_eval_model(tX,P);
            Q(j) = plot((max_scaleT/maxT).*tX,X,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s: %.1f min, %.1f hr, %.1f hr, %.1f]',...
                name_r{j},gc{k},60*log(2)/P(2),P([3 4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot((max_scaleT/maxT).*u,s,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s]',name_r{j},gc{k});
        end
        %z = D(k,iP);
        %plot(x./mxS,z,'s','markersize',10,'color',plotcolor);
        %[u,~,w] = unique(x);
        %s = accumarray(w,z)./accumarray(w,1);
        %plot(u./mxS,s,'--','linewidth',2,'color',plotcolor);
    end
end
