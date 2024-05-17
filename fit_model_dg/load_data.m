function [D,t,iM,iP,gid,tid,gc] = load_data(name)

X = importdata(name);
gid = X.textdata(2:end,3);
tid = X.textdata(2:end,2);
gc = X.textdata(2:end,1);
D = X.data;

t = X.textdata(1,4:end);
iM = strncmp(t,'transcript_',11) + strncmp(t,'M_',2) > 0;
iP = strncmp(t,'precursor_',10) + strncmp(t,'P_',2) > 0;
t = regexprep(t,'transcript_','');
t = regexprep(t,'M_','');
t = regexprep(t,'precursor_','');
t = regexprep(t,'P_','');
t = regexprep(t,'\"','');
t = regexprep(t,'_.*','');
t = str2double(t);
