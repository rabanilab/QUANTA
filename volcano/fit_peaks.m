function [X,Xu,Xw] = fit_peaks(Sid,S,Kid,Kweight,output_pref)
% Output:
%  X  = all peaks: <peak weight> <reporter id> <start> <end> <sequence>
%       peak weight = sum all weights of peak positions / peak size
%  Xu = unique peaks: <peak weight> <sequence>
%  Xw = unique peak sequences: <peak weight> <sequence>

if (nargin < 5)
    output_pref = [];
end

% positional weights
weight_prctile = 99;

% peak selection
window_size = 8;
nan_pct = 0.2;
peak_select_prctile = 75;
min_peak_width = 3;
peak_extend = 1;


% ------------------------
% positional weights
% ------------------------
W = cell(size(S));
for j = 1:size(S,1)
    W{j} = assign_positional_weights(S{j},Kid,Kweight);
end
w = convert_weights_matrix(W);
pmin = prctile(w(w>0),100-weight_prctile);
pmax = prctile(w(w>0),weight_prctile);
if (pmax-pmin < 1)
    pmin = min(w(:));
    pmax = max(w(:));
end
wnn = (w-pmin)./(pmax-pmin);
wnn(wnn>1) = 1;
wnn(wnn<0) = 0;
fprintf('[1] Positional weights: pct = %.1f [%.2f %.2f]\n',weight_prctile,pmin,pmax);


% ------------------------
% extract peak sequences
% ------------------------

% Calculate a threshold for peak selection
% use a sliding window average percentile
d2 = [];
maxN = find(sum(~isnan(wnn),1)./size(wnn,1)>=nan_pct,1,'last');
for j = 1:(maxN-window_size+1)
    d2 = [d2 nanmean(wnn(:,j:(j+window_size-1)),2)];
end
t1 = prctile(d2(d2>0),peak_select_prctile);
fprintf('[1] Peak selection: pct = %.1f, thr = %.2f\n',peak_select_prctile,t1);

% select peaks
[X,Xu,Xw] = peaks_extract(Sid,S,wnn,t1,min_peak_width,peak_extend);
if (~isempty(output_pref))
    write_text_file([output_pref '.a.txt'],X);
    write_text_file([output_pref '.u.txt'],Xu);
    write_text_file([output_pref '.w.txt'],Xw);
end



% --------------------------------------------------------------------
% kmer positional weights
% --------------------------------------------------------------------

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
sc = log2(sc+1);

function W = convert_weights_matrix(X)

Wlen = cellfun(@length,X);
W = nan(size(X,1),max(Wlen));
for i = 1:size(X,1)
    W(i,1:length(X{i})) = X{i};
end


% --------------------------------------------------------------------
% identify peaks
% --------------------------------------------------------------------

function [X,Xu,Xw] = peaks_extract(id,seq,w,min_peak_weight,min_peak_width,peak_extend)
% find peaks in input weights
%
% Output:
%  X  = all peaks: <peak weight> <reporter id> <start> <end> <sequence>
%       peak weight = sum all weights of peak positions / peak size
%  Xu = unique peaks: <peak weight> <sequence>
%  Xw = unique peak sequences: <peak weight> <sequence>

w(isnan(w)) = 0;
rid = repmat((1:size(w,1))',1,size(w,2))';
cid = repmat((1:size(w,2)),size(w,1),1)';
q = find(w' >= min_peak_weight)';

X = [];
Xu = [];
Xw = [];
if (isempty(q))
    return;
end

% find peaks in input weights
s = q([0 q(2:end) == q(1:end-1)+1]==0)';
e = q([q(2:end) == q(1:end-1)+1 0]==0)';
if (sum(rid(s)~=rid(e))==0)
    
    % extend peaks
    P = [rid(s) cid(s)-peak_extend cid(e)+peak_extend];
    
    % extract peak sequences and weights
    m = size(P,1);
    S = cell(m,1);
    W = zeros(m,1);
    for k = 1:m
        if (P(k,2) < 1)
            P(k,2) = 1;
        end
        if (P(k,3) > length(seq{P(k,1)}))
            P(k,3) = length(seq{P(k,1)});
        end
        S{k} = seq{P(k,1)}(P(k,2):P(k,3));
        W(k) = sum(w(P(k,1),P(k,2):P(k,3)))./(P(k,3)-P(k,2)+1);
    end
    X = [num2cell(W) id(P(:,1)) num2cell(P(:,2:end)) S];

    % filter by peak length
    l = cellfun(@length,S);
    j = (l >= min_peak_width+2*peak_extend);
    X = X(j==1,:);
end

np = size(X,1);
fprintf('Peaks: %d peaks identified\n', np);

if (np > 0)
    % filter peaks with same score and same sequence
    Xu = filter_unique_peaks(X(:,[1 5]));

    % merge peaks with the same sequence and sum their scores
    %Xw = merge_peaks_sequences(Xu);
    Xw = merge_peaks_sequences(X(:,[1 5]));

    plen = cellfun(@length,Xw(:,2));
    pscore = abs(cell2mat(Xw(:,1)));
    fprintf('Peaks: %d unique peaks (sequence+score) identified (%.1f%%)\n',size(Xu,1),100*size(Xu,1)./np);
    fprintf('Peaks: %d unique peak sequences identified (%.1f%%)\n', size(Xw,1),100*size(Xw,1)./np);
    fprintf('Peaks: mean length = %.1f (min=%d, max=%d)\n', mean(plen),min(plen),max(plen));
    fprintf('Peaks: mean absolute height = %.2f (min=%.2f, max=%.2f)\n', mean(pscore),min(pscore),max(pscore));
    fprintf('Peaks: %.1f%% of bases are under peaks\n', 100*sum(cellfun(@length,X(:,5)))./numel([seq{:}]));
else
    Xu = [];
    Xw = [];
end

function U = filter_unique_peaks(Z)
% filter peaks with same score and same sequence
% U = [score] [sequence]

W(:,1) = cell2mat(Z(:,1));
[u,~,t] = unique(Z(:,2));
W(:,2) = t;
W = unique(W,'rows');

U = [num2cell(W(:,1)) u(W(:,2))];
U = sortrows(U,-1);

function U = merge_peaks_sequences(Z)
% merge peaks with the same sequence and sum their scores
% U = [score] [sequence]

[W1,~,t] = unique(Z(:,2));
W2 = accumarray(t,cell2mat(Z(:,1)));
U = [num2cell(W2) W1];
U = sortrows(U,-1);

