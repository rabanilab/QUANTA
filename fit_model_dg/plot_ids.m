function plot_ids(all_ids,all_data_a,all_data_r,pref,run_mfit,gene_ids,mzt,feps,prefix)

if (nargin < 5)
    run_mfit = 1;
end
if (nargin < 6)
    gene_ids = 0;
end
if (nargin < 7)
    mzt = [];
end
if (nargin < 8)
    feps = 0;
end
if (nargin < 9)
    prefix = 'data_from_drive/workdir/exp_df_classified_';
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
for i = 1:max(size(all_ids))
    clf;
    C = colormap('lines');
    
    % polyA datasets
    subplot(1,2,1);
    hold on;
    [gname,Q,S,mxT] = plot_a(all_ids{i},all_data_a,prefix,run_mfit,gene_ids,mzt,C);
    hold off;
    xlabel('time (hr)');
    ylabel('logFPKM');
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([gname ' (polyA)']);
    set(gca,'ylim',[-3 15],'xlim',[0 mxT+0.5],'fontsize',15);
    if (feps > 0)
        set(gca,'ylim',[-3 10]);
        axis square;
        legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    end
    
    % ribo depletion datasets
    subplot(1,2,2);
    hold on;
    [gname,Q,S,mxT] = plot_r(all_ids{i},all_data_r,prefix,run_mfit,gene_ids,mzt,C);
    hold off;
    xlabel('time (hr)');
    ylabel('logFPKM');
    legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    title([gname ' (ribo)']);
    set(gca,'ylim',[-3 15],'xlim',[0 mxT+0.5],'fontsize',15);
    if (feps > 0)
        set(gca,'ylim',[-3 10]);
        axis square;
        legend(Q(Q>0),S(Q>0),'box','off','location','northoutside');
    end
    
    % save file
    if (~isempty(gname))
        saveas(h, [pref gname '.jpg'],'jpg');
        if (feps > 0)
            saveas(h, [pref gname '.eps'],'epsc');
        end
    else
        saveas(h, [pref all_ids{i} '.EMPTY.jpg'],'jpg');
    end
end

close all;


function [gname,Q,S,mxT] = plot_a(id,all_data_a,prefix,run_mfit,gene_ids,mzt,C)

name_a = regexprep(all_data_a,'_',' ');

S = {};
Q = [];
gname = [];
mxT = 1;
for j = 1:max(size(all_data_a))
    [D,t,iM,iP,gid,tid,gc] = load_data([prefix all_data_a{j} '.txt']);
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
        plot(x,y,'.','markersize',35,'color',C(j,:));
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
            Q(j) = plot(tX,X,'-','linewidth',2,'color',C(j,:));
            S{j} = sprintf('%s [%s: %.1f min, %.1f min, %.1f hr, %.1f]',...
                name_a{j},gc{k},60*log(2)/P(2),60*log(2)/P(3),P([4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot(u,s,'-','linewidth',2,'color',C(j,:));
            S{j} = sprintf('%s [%s]',name_a{j},gc{k});
        end
        z = D(k,iP);
        plot(x,z,'s','markersize',10,'color',C(j,:));
        [u,~,w] = unique(x);
        s = accumarray(w,z)./accumarray(w,1);
        plot(u,s,'--','linewidth',2,'color',C(j,:));
    end
end


function [gname,Q,S,mxT] = plot_r(id,all_data_r,prefix,run_mfit,gene_ids,mzt,C)

name_r = regexprep(all_data_r,'_',' ');

S = {};
Q = [];
gname = [];
mxT = 1;
for j = 1:max(size(all_data_r))
    [D,t,iM,iP,gid,tid,gc] = load_data([prefix all_data_r{j} '.txt']);
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
        plot(x,y,'.','markersize',35,'color',C(j,:));
        if (run_mfit)
            if (size(x,2)>=4)
                P = fit_model(x,y,mzt(2:3));
            else
                P = [mean(y) 0 0 0];
            end
            X = dg_eval_model(tX,P);
            Q(j) = plot(tX,X,'-','linewidth',2,'color',C(j,:));
            S{j} = sprintf('%s [%s: %.1f min, %.1f hr, %.1f hr, %.1f]',...
                name_r{j},gc{k},60*log(2)/P(2),P([3 4 1]));
        else
            [u,~,w] = unique(x);
            s = accumarray(w,y)./accumarray(w,1);
            Q(j) = plot(u,s,'-','linewidth',2,'color',C(j,:));
            S{j} = sprintf('%s [%s]',name_r{j},gc{k});
        end
        z = D(k,iP);
        plot(x,z,'s','markersize',10,'color',C(j,:));
        [u,~,w] = unique(x);
        s = accumarray(w,z)./accumarray(w,1);
        plot(u,s,'--','linewidth',2,'color',C(j,:));
    end
end
