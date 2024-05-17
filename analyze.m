% ---------------------------------------------------------------------
% preparation
% ---------------------------------------------------------------------

% load code
addpath('fit_model_dg');
addpath('volcano');
addpath('misc');

% data
all_data_r = {'Medina2021_ribo'}
all_data_a = {'Medina2021_polya'}
all_data_c = all_data_r;
mzt_class = 3.5;
mzt_fit = [2 6 6]; % [min dg onset time; 2p model only] [max dg onset time] [min dg offset time]
kmer_dir = 'volcano/gene_kmers';
seq_file = 'volcano/gene_kmers/3utr_seqs.genes.txt';
data_dir = 'data_files';

% prepare results dir
mkdir('results');

% ---------------------------------------------------------------------
% create kmer counts data
% ---------------------------------------------------------------------
kmer_run = 1;

if (kmer_run)
    parfor k = 1:7
        extract_kmer_counts(k,seq_file,[kmer_dir '/kmers_counts.' num2str(k)]);
    end
end

% ---------------------------------------------------------------------
% create parameter XML file
% ---------------------------------------------------------------------
xml_create = 1;

if (xml_create)
    create_param_xml(all_data_r,all_data_a,all_data_c,mzt_class,mzt_fit,...
    kmer_dir,seq_file,data_dir,0.75,[1 0 0 0]);
end

% ---------------------------------------------------------------------
% Modeling analysis
% ---------------------------------------------------------------------
model_analysis = 1;

if (model_analysis)
    analyze_temporal_data();
end
