function extract_kmer_counts(W,sequence_file,filename)

fprintf('Extracting: file=%s\n',sequence_file);

f = fopen(sequence_file);
X = textscan(f,'%s %s');
fclose(f);
id = X{1};
seq = char(upper(X{2}));
seq = cellstr(seq);
fprintf('Sequences: n=%d\n',size(seq,1));

% count k-mers
for K = W
    fprintf('Extracting: k=%d\n',K);
    
    Xrows = all_kmers(K);
    Xcols = id;
    X = sparse(size(Xrows,1),size(Xcols,1));
    for i = 1:size(seq,1)
        ci = nmercount(seq{i},K,0);
        if (isempty(Xrows))
            Xrows = ci(:,1);
            X = [X cell2mat(ci(:,2))];
        else
            fprintf('Counts: k=%d, n=%d\n',K,i);
            [x,y] = ismember(Xrows,ci(:,1)); % [Xrows(x) ci(y(y>0),1)]
            X(x,i) = cell2mat(ci(y(y>0),2));
        end
    end
    
    % remove empty k-mers
    i = sum(X,2) > 0;
    fprintf('size %d K-mers (%d sequences), %d K-mers exist\n', K, size(id,1), full(sum(i)));
    X = X(i,:);
    Xrows = Xrows(i);
    
    save([filename '.' num2str(K) '.mat'],'X','Xrows','Xcols');
end

function S = all_kmers(K)

L = ['A'; 'C'; 'G'; 'T'];

S = [];
if (K == 1)
    S = cellstr(L);
elseif (K > 1)
    Si = all_kmers(K-1);
    for i = 1:max(size(L))
        S = [S;strcat(L(i),Si)];
    end
end
