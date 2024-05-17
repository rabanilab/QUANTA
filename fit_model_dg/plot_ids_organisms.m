function plot_ids_organisms(ids,data_a1,data_r1,loadpref1,mzt1,outpref,scale_time,run_mfit,gene_ids,feps)

if (nargin < 7)
    scale_time = 0;
end
if (nargin < 8)
    run_mfit = 1;
end
if (nargin < 9)
    gene_ids = 0;
end
if (nargin < 10)
    feps = 0;
end

N = max(size(data_a1));

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:max(size(ids))
    clf;
    C = colormap('lines');

    % polyA datasets
    subplot(1,2,1);
    hold on;
    mxT = 0;
    Q = [];
    S = [];
    T = [];
    for j = 1:N
        [gname1,Q1,S1,mxT1] = plot_a(ids{i,j},data_a1{j},loadpref1{j},scale_time,run_mfit,gene_ids,mzt1{j},C(j,:));
        mxT = max(mxT,mxT1);
        Q = [Q Q1];
        S = [S S1];
        if (isempty(T))
            T = gname1;
        else
            T = [T ',' gname1];
        end
    end
    hold off;
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([T ' (polyA)']);
    ylabel('logFPKM');
    if (scale_time)
        xlabel('scaled time (a.u.)');
        set(gca,'ylim',[-3 15],'xlim',[0 1],'fontsize',15);
    else
        xlabel('time (hr)');
        set(gca,'ylim',[-3 15],'xlim',[0 mxT+0.5],'fontsize',15);
    end
    if (feps > 0)
        %set(gca,'ylim',[-3 15]);
        axis square;
        legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    end

    % ribo depletion datasets
    subplot(1,2,2);
    hold on;
    mxT = 0;
    Q = [];
    S = [];
    T = [];
    for j = 1:N
        [gname1,Q1,S1,mxT1] = plot_r(ids{i,j},data_r1{j},loadpref1{j},scale_time,run_mfit,gene_ids,mzt1{j},C(j,:));
        mxT = max(mxT,mxT1);
        Q = [Q Q1];
        S = [S S1];
        if (isempty(T))
            T = gname1;
        else
            T = [T ',' gname1];
        end
    end
    hold off;
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([T ' (ribo)']);
    ylabel('logFPKM');
    if (scale_time)
        xlabel('scaled time (a.u.)');
        set(gca,'ylim',[-3 15],'xlim',[0 1],'fontsize',15);
    else
        xlabel('time (hr)');
        set(gca,'ylim',[-3 15],'xlim',[0 mxT+0.5],'fontsize',15);
    end
    if (feps > 0)
        %set(gca,'ylim',[-3 15]);
        axis square;
        legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    end

    % save file
    if (~isempty(gname1))
        saveas(h, [outpref gname1 '.jpg'],'jpg');
        if (feps > 0)
            saveas(h, [outpref gname1 '.eps'],'epsc');
        end
    else
        saveas(h, [outpref ids{i} '.EMPTY.jpg'],'jpg');
    end
end

close all;


function [gname,Q,S,mxT] = plot_a(id,data_a,loadpref,scale_time,run_mfit,gene_ids,mzt,plotcolor)

name_a = regexprep(data_a,'_',' ');

mxS = 1;
if (scale_time)
    for j = 1:max(size(data_a))
        [~,t] = load_data([loadpref data_a{j} '.txt']);
        if (max(t)>mxS)
            mxS = max(t);
        end
    end
end

S = {};
Q = [];
gname = [];
mxT = 10;
for j = 1:max(size(data_a))
    [D,t,iM,iP,gid,tid,gc] = load_data([loadpref data_a{j} '.txt']);
    if (max(t)>mxT)
        mxT = max(t);
    end
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
        plot(x./mxS,y,'.','markersize',35,'color',plotcolor);
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
            Q(j) = plot(tX./mxS,X,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s: %.1f min, %.1f min, %.1f hr, %.1f]',...
                name_a{j},gc{k},60*log(2)/P(2),60*log(2)/P(3),P([4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot(u./mxS,s,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s]',name_a{j},gc{k});
        end
        %z = D(k,iP);
        %plot(x./mxS,z,'s','markersize',10,'color',plotcolor);
        %[u,~,w] = unique(x);
        %s = accumarray(w,z)./accumarray(w,1);
        %plot(u./mxS,s,'--','linewidth',2,'color',plotcolor);
    end
end


function [gname,Q,S,mxT] = plot_r(id,data_r,loadpref,scale_time,run_mfit,gene_ids,mzt,plotcolor)

name_r = regexprep(data_r,'_',' ');

mxS = 1;
if (scale_time)
    for j = 1:max(size(data_r))
        [~,t] = load_data([loadpref data_r{j} '.txt']);
        if (max(t)>mxS)
            mxS = max(t);
        end
    end
end

S = {};
Q = [];
gname = [];
mxT = 10;
for j = 1:max(size(data_r))
    [D,t,iM,iP,gid,tid,gc] = load_data([loadpref data_r{j} '.txt']);
    if (max(t)>mxT)
        mxT = max(t);
    end
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
        plot(x./mxS,y,'.','markersize',35,'color',plotcolor);
        if (run_mfit)
            if (size(x,2)>=4)
                P = fit_model(x,y,mzt(2:3));
            else
                P = [mean(y) 0 0 0];
            end
            X = dg_eval_model(tX,P);
            Q(j) = plot(tX./mxS,X,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s: %.1f min, %.1f hr, %.1f hr, %.1f]',...
                name_r{j},gc{k},60*log(2)/P(2),P([3 4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot(u./mxS,s,'-','linewidth',2,'color',plotcolor);
            S{j} = sprintf('%s [%s]',name_r{j},gc{k});
        end
        %z = D(k,iP);
        %plot(x./mxS,z,'s','markersize',10,'color',plotcolor);
        %[u,~,w] = unique(x);
        %s = accumarray(w,z)./accumarray(w,1);
        %plot(u./mxS,s,'--','linewidth',2,'color',plotcolor);
    end
end
