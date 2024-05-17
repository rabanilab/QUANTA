function [mT,mD,mM,mM0] = dg_eval_all(all_data,mtid,mP,model,datafile,pref)

n = max(size(all_data));

mT = cell(1,n); % timepoints
mD = cell(1,n); % data
mM = cell(1,n); % model predictions
mM0 = cell(1,n); % logX0 values (null model ?)
for j = 1:n
    if (isempty(datafile))
        [D,t,iM,~,~,tid] = load_data([pref all_data{j} '.txt']);
    else
        [D,t,iM,~,~,tid] = load_data(datafile);
    end
    [k1,k2] = ismember(mtid,tid); %[tid(k2(k2>0) mtid(k1))]
    x = t(iM);
    y = D(k2(k2>0),iM);
    mD{j} = nan(size(mtid,1),size(y,2));
    mD{j}(k1,:) = y;
    mT{j} = x;
    mM0{j} = repmat(mP{j}(:,1),1,size(x,2));
    for i = 1:size(mP{j},1)
        if (model == 2)
            mM{j}(i,:) = dg_eval_model_2p(x, mP{j}(i,:));
        else
            mM{j}(i,:) = dg_eval_model(x, mP{j}(i,:));
        end
    end
end