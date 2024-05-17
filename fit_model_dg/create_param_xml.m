function create_param_xml(all_data_r,all_data_a,all_data_c,mzt_class,mzt_fit,...
    kmer_dir,seq_file,data_dir, ...
    class_FC,organism_array,is_mzt,cnt_gene_norm,...
    fpkm_pref,gene_file,list_file)
% class_FC = fold-reduction after/before that is smaller than this value is excluded
% set to 0.5 for polyA/non-MZT datasets
% set to 0.75 for total-RNA/MZT datasets
% organism_array = [run_zfish run_frog run_mouse run_human]

if (nargin < 9)
    class_FC = 0.75;
end
if (nargin < 10)
    organism_array = [0 0 0 0];
end
if (nargin < 11)
    is_mzt = 1;
end
if (nargin < 12)
    cnt_gene_norm = [1 0];
end
if (nargin < 13)
    fpkm_pref = [data_dir '/exp_df_classified_'];
end
if (nargin < 14)
    gene_file = [data_dir '/mart_export_filterd.txt'];
end
if (nargin < 15)
    list_file = [data_dir '/file_list.txt'];
end

organism_ids = {'run_zfish' 'run_frog' 'run_mouse' 'run_human'};

% open xml doc
docNode = com.mathworks.xml.XMLUtils.createDocument('param');
entry_node = docNode.createElement('Entry');

% write param
node = docNode.createElement('kmer_dir');
nname = docNode.createTextNode(kmer_dir);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
node = docNode.createElement('seq_file');
nname = docNode.createTextNode(seq_file);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
node = docNode.createElement('data_dir');
nname = docNode.createTextNode(data_dir);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
node = docNode.createElement('fpkm_pref');
nname = docNode.createTextNode(fpkm_pref);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
node = docNode.createElement('gene_file');
nname = docNode.createTextNode(gene_file);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
node = docNode.createElement('list_file');
nname = docNode.createTextNode(list_file);
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);

node = docNode.createElement('class_FC');
nname = docNode.createTextNode(num2str(class_FC));
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
for i = 1:4
    node = docNode.createElement(organism_ids{i});
    nname = docNode.createTextNode(num2str(organism_array(i)));
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end
node = docNode.createElement('is_mzt');
nname = docNode.createTextNode(num2str(is_mzt));
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
for i = 1:max(size(cnt_gene_norm))
    node = docNode.createElement(['cnt_gene_norm' num2str(i)]);
    nname = docNode.createTextNode(num2str(cnt_gene_norm(i)));
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end

node = docNode.createElement('mzt_class');
nname = docNode.createTextNode(num2str(mzt_class));
node.appendChild(nname);
docNode.getDocumentElement.appendChild(node);
for i = 1:max(size(mzt_fit))
    node = docNode.createElement('mzt_fit');
    nname = docNode.createTextNode(num2str(mzt_fit(i)));
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end
for i = 1:max(size(all_data_r))
    node = docNode.createElement('data_r');
    nname = docNode.createTextNode(all_data_r{i});
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end
for i = 1:max(size(all_data_a))
    node = docNode.createElement('data_a');
    nname = docNode.createTextNode(all_data_a{i});
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end
for i = 1:max(size(all_data_c))
    node = docNode.createElement('data_c');
    nname = docNode.createTextNode(all_data_c{i});
    node.appendChild(nname);
    docNode.getDocumentElement.appendChild(node);
end

% create file
text = xmlwrite(docNode);
fid = fopen('param.xml','wt');
fprintf(fid, '%s\n', text);
fclose(fid);
