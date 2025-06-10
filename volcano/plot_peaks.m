function h = plot_peaks(xid,Sid,S,p,kmer_pv,L)
% p=2 (negative kmers),p=3 (positive kmers)

j = find(strcmp(xid,Sid));
if (isempty(j))
    h = [];
    return;
end

n = max(size(kmer_pv));
W = [];
for i = 1:n
    Kid = kmer_pv{i}(:,1);
    Kweight = -1*log10(cell2mat(kmer_pv{i}(:,p)));    
    wi = assign_positional_weights(S{j},Kid,Kweight);
    W{i} = wi;
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
mx = -inf;
hold on;
for i = 1:n
    l = length(W{i});
    plot(1:l,W{i}(1:l)','LineWidth',3);
    mx = max(mx,max(W{i}(1:l)));
end
hold off;
axis tight;
set(gca,'xtick',1:l,'xticklabel',cellstr(S{j}')');
mx = max(mx,50);
set(gca,'ylim',[0 mx],'fontsize',18);
title(xid);
ylabel('weight(log2)');
legend(L,'box','off');

function sc = assign_positional_weights(s,Kstr,Kpval)

K_len = cellfun(@length,Kstr);
sc = zeros(size(s));
for k = unique(K_len)'
    for x = find(K_len==k)'
        p = regexp(s,Kstr{x}); % all positions within s
        for y = 0:k-1
            q = p+y;
            q = q(q<size(s,2));
            sc(q) = sc(q) + 10^Kpval(x);
        end
    end
end
sc = log2(sc+1) - 10;
sc(sc<0) = 0;

