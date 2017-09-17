function linearmodels = train_linearmodel_in_leaves(model, xTrain_c_all, yTrain_all, numb,nxc)

xTrain_d = model.X;
NTrees = model.NTrees;

h0 = waitbar(0, 'Synthesis Forest');

for i = 1:NTrees
    
    leaf_index = find((model.Trees{i}.Children(:,1)==0)&(model.Trees{i}.Children(:,2)==0));
    numleafs = length(leaf_index);
    [~,node] = predict(model.Trees{i}, xTrain_d(~model.OOBIndices(:,i),:));
    xTrain_c = xTrain_c_all(~model.OOBIndices(:,i),:);
    yTrain = yTrain_all(~model.OOBIndices(:,i),:);
    h = waitbar(0, 'Sythesis Tree');
    
    for j=1:numleafs
        
        % find indices of samples which end up in this leaf
        leafmodels(j).leaves = {find(node==leaf_index(j))}; %#ok<*AGROW>
        
        % mean prediction at this leaf
        leafmodels(j).mean = model.Trees{i}.NodeMean(leaf_index(j));
        
        % the control variables sample values which contribute to this leaf (support)
        leafmodels(j).xdata = {xTrain_c(leafmodels(j).leaves{1,1},:)};
        
        % the response variable value which contribute to this leaf
        leafmodels(j).ydata = {yTrain(leafmodels(j).leaves{1,1})};
        
        % train a linear model
        X = [ones(size(leafmodels(j).ydata{1},1),1), leafmodels(j).xdata{1}];
        Y = leafmodels(j).ydata{1,1};
        
        Y_t = Y;
        X_t = X(:,1:nxc*numb+1);
        if size(X_t,1)<size(X_t,2)
            % keyboard;
            leafmodels(j).coeff = nan(nxc*numb+1,1);
        else
            options=optimset;
            options.Algorithm='interior-point'; % to get rid of some warnings
            options.Display='off';
            Aineq = diag([0;-ones(nxc*numb,1)]);
            bineq = zeros(nxc*numb+1,1);
            lb=zeros(1,nxc*numb+1);
            coeff = lsqlin(X_t,Y_t,Aineq,bineq,[],[],lb,[],[],options);
            leafmodels(j).coeff = coeff;
            leafmodels(j).xdata = [];
            leafmodels(j).ydata = [];
            leafmodels(j).leaves = [];
        end
        progress = j/numleafs;
        waitbar(progress, h, sprintf('Leaf %d of %d', j, numleafs));
    end
    
    linearmodels{i} = leafmodels;
    close(h);
    
    progress = i/NTrees;
    waitbar(progress, h0, sprintf('Tree %d of %d', i, NTrees));
end

close(h0)