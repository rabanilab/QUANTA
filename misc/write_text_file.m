function write_text_file(filename, data)
% write text data into a file
% input is a cell

rows = size(data,1);
cdata = to_str(data);

fid = fopen(filename, 'wt');
for (i=1:rows)
    fprintf(fid, '%s\t', cdata{i,1:end-1});
    fprintf(fid, '%s\n', cdata{i,end});
end
fclose(fid);



function h = to_str(v)

[n,m] = size(v);
h = cell(n,m);
for i=1:n
    for j=1:m
        if(isnumeric(v{i,j})+islogical(v{i,j})>0)
            h{i,j} = num2str(v{i,j});
        elseif (isempty(v{i,j}))
            h{i,j} = '';
        else
            h{i,j} = v{i,j};
        end
    end
end
