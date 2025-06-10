function h = plot_volcano(e,p1,p2,kmer_ids,alphaT,effectT,text_lim,e_lim,p_lim)

if (nargin < 5)
    alphaT = 0.01;
end
if (nargin < 6)
    effectT = prctile(abs(e),80);
end
if (nargin < 7)
    text_lim = 20;
end
if (nargin < 8)
    e_lim = [-1 1];
end
if (nargin < 9)
    p_lim = [0 30];
end

NORM_VAR = 1e-5;
if (max(size(unique(e(:))))<=1)
    e = e + NORM_VAR*randn(size(e));
end
if (max(size(unique(p1(:))))<=1)
    p1 = p1 + NORM_VAR*randn(size(p1));
end
if (max(size(unique(p2(:))))<=1)
    p2 = p2 + NORM_VAR*randn(size(p2));
end

% effect size
x1 = e';
x1(x1<e_lim(1)) = e_lim(1);
x1(x1>e_lim(2)) = e_lim(2);

% pvalue
x2 = -1*log10(nanmin(p1,p2)');
x2(x2<p_lim(1)) = p_lim(1);
x2(x2>p_lim(2)) = p_lim(2);

% add random noise for plotting purposes
k = (~isnan(x2)).*(~isnan(x1)) == 1;
if (max(size(unique(x1(k))))<2)
    x1(k) = x1(k) + randn(sum(k),1)*0.01;
end
if (max(size(unique(x2(k))))<2)
    x2(k) = x2(k) + randn(sum(k),1)*0.01;
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) 0.5*scrsz(3) scrsz(4)]);

generate_plot(x1,x2,kmer_ids,-1*log10(alphaT),effectT,text_lim);
%generate_plot2(x1,x2,kmer_ids,-1*log10(alphaT),effectT,text_lim,e_lim,p_lim);


function [k0,k1] = generate_plot(E,P,kmer_ids,alphaT,effectT,text_lim)

k0 = isnan(P) + isnan(E) == 0;
k1 = (P >= alphaT).*(abs(E) >= effectT) == 1;
mxPj = max(max(P),alphaT) + 1;
mxEj = max(max(abs(E)),effectT) + 0.1;

clf;
hold on;
if ((sum(k0)>1))
    s = scatter(E(k0),P(k0),'k','MarkerFaceColor','k');
    alpha(s,0.01);

    if (text_lim > 0)
        write_text(P,E,kmer_ids,text_lim,alphaT,effectT);
    end
    if (sum(k1) > 1)
        n1 = zeros(size(E(k1)));
        n2 = zeros(size(P(k1)));
        if (max(size(unique(E(k1))))<2)
            n1 = 0.01*randn(size(E(k1)));
        end
        if (max(size(unique(P(k1))))<2)
            n2 = 0.01*randn(size(P(k1)));
        end
        dscatter(E(k1)+n1,P(k1)+n2,'MSIZE',50,'MARKER','o');
    else
        plot(E(k1),P(k1),'.','markerSize',20);
    end
    line([0 0],[0 mxPj],'LineStyle','-','color','k','linewidth',2);
    line([-1*mxEj mxEj],[alphaT alphaT],'LineStyle','-','color','r','linewidth',1);
    line([effectT effectT],[0 mxPj],'LineStyle','-','color','r','linewidth',1);
    line(-1*[effectT effectT],[0 mxPj],'LineStyle','-','color','r','linewidth',1);
    axis tight;
    axis square;
    box on;
    set(gca,'xlim',[-1*mxEj mxEj],'ylim',[0 mxPj],'fontsize',12);
end
hold off;
xlabel('effect size (SMD)');
ylabel('p-value');
title(sprintf('n=%d (all), n=%d (%.1f)',sum(k0),sum(k1),100*sum(k1)/sum(k0)));


function generate_plot2(E,P,kmer_ids,alphaT,effectT,text_lim,e_lim,p_lim)

clf;
hold on;
k = (~isnan(P)).*(~isnan(E)) == 1;
if ((sum(k)>1))
    %s = plot(x1(k),x2(k),'.k','markerSize',20);
    %alpha(s,.5);
    dscatter(E(k),P(k),'MSIZE',50,'MARKER','o');
    line(e_lim,[alphaT alphaT],'LineStyle','-','color','k','linewidth',1);
    line(-1*[effectT effectT],p_lim,'LineStyle','-','color','k','linewidth',1);
    line([effectT effectT],p_lim,'LineStyle','-','color','k','linewidth',1);
    if (text_lim > 0)
        write_text(P,E,kmer_ids,text_lim,alphaT,effectT);
    end
end
hold off;
xlabel('effect size');
ylabel('-log(p-value)');
set(gca,'xlim',e_lim,'ylim',p_lim,'fontsize',12);


function write_text(P,E,kmer_ids,text_lim,alphaT,effectT)

k = (P >= alphaT).*(abs(E) >= effectT) == 1;
s2 = sort(P(k));
if (~isempty(s2))
    lastp = text_lim;
    if (size(s2,1)<text_lim)
        lastp = size(s2,1);
    end
    alphaT_new = max(s2(end-lastp+1),alphaT);

    q = (P >= alphaT_new).*(abs(E) >= effectT) == 1;
    s = find(q);
    if (max(size(s))>text_lim)
        s = s(1:text_lim);
    end
    for i = s'
        seq = regexprep(kmer_ids(i),'T','U');
        text(E(i),P(i),seq,'fontsize',10);
    end
end

