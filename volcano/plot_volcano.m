function h = plot_volcano(e,p1,p2,kmer_ids,alphaT,effectT,maxText,e_lim,p_lim)

if (nargin < 5)
    alphaT = 0.01;
end
if (nargin < 6)
    effectT = prctile(abs(e),80);
end
if (nargin < 7)
    maxText = 20;
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

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) 0.5*scrsz(3) scrsz(4)]);

x1 = e';
x1(x1<e_lim(1)) = e_lim(1);
x1(x1>e_lim(2)) = e_lim(2);
x2 = -1*log(nanmin(p1,p2)')./log(10);
x2(x2<p_lim(1)) = p_lim(1);
x2(x2>p_lim(2)) = p_lim(2);
x0 = -1*log(alphaT)./log(10);
k = (~isnan(x2)).*(~isnan(x1)) == 1;
if ((sum(k)>1))
    if (max(size(unique(x1(k))))<2)
        x1(k) = x1(k) + randn(sum(k),1)*0.01;
    end
    if (max(size(unique(x2(k))))<2)
        x2(k) = x2(k) + randn(sum(k),1)*0.01;
    end
    hold on;
    dscatter(x1(k),x2(k),'MSIZE',50,'MARKER','o');
    line(e_lim,[x0 x0],'LineStyle','-','color','k','linewidth',1);
    line(-1*[effectT effectT],p_lim,'LineStyle','-','color','k','linewidth',1);
    line([effectT effectT],p_lim,'LineStyle','-','color','k','linewidth',1);
    hold off;
    
    s2 = sort(x2(k.*(abs(x1)>effectT)==1));
    if (~isempty(s2))
        lastp = maxText;
        if (size(s2,1)<maxText)
            lastp = size(s2,1);
        end
        d = max(s2(end-lastp+1),x0);
        s = find(k.*(x2>=d).*(abs(x1)>effectT)==1);
        if (max(size(s))>maxText)
            s = s(1:maxText);
        end
        for i = s'
            text(x1(s),x2(s),kmer_ids(s),'fontsize',14);
        end
    end
end
xlabel('effect size');
ylabel('-log(p-value)');
set(gca,'xlim',e_lim,'ylim',p_lim,'fontsize',14);

