function [W1,W2,Ethr] = volcano_mw(gid,rates,kname,Krange,ESIZE,ALPHA,Kdir,e_lim,p_lim,MHcorrect)

if (nargin < 4)
    Krange = 3:7;
end
if (nargin < 5)
    ESIZE = 0.3;
end
if (nargin < 6)
    ALPHA = 0.01;
end
if (nargin < 7)
    Kdir = '/Users/mrabani/GoogleDrive/Matlab/MATLAB/volcano/utrseq_kmers';
end
if (nargin < 8)
    e_lim = [-1 1];
end
if (nargin < 9)
    p_lim = [0 30];
end
if (nargin < 10)
    MHcorrect = 'fdr';
end

[W1,W2,Ethr] = volcano_test(gid,rates,kname,Krange,ESIZE,ALPHA,Kdir,e_lim,p_lim,'mw',MHcorrect);


