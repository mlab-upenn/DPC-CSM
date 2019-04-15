function linearmodels = train_linearmodel_in_leaves(model, xTrain_c, yTrain, idm)

leaf_index = find((model.Children(:,1)==0)&(model.Children(:,2)==0));
numleafs = length(leaf_index);
[~,node] = resubPredict(model);

h = waitbar(0, 'Sythesis Tree');

for idl=1:numleafs

    % find indices of samples which end up in this leaf
    linearmodels(idl).leaves = {find(node==leaf_index(idl))}; %#ok<*AGROW>

    % mean prediction at this leaf
    linearmodels(idl).mean = model.NodeMean(leaf_index(idl));

    % the control variables sample values which contribute to this leaf (support)
    linearmodels(idl).xdata = {xTrain_c(linearmodels(idl).leaves{1,1},:)};

    % the response variable value which contribute to this leaf
    linearmodels(idl).ydata = {yTrain(linearmodels(idl).leaves{1,1})};

    % train a linear model
    X = [ones(size(linearmodels(idl).ydata{1},1),1), linearmodels(idl).xdata{1}];
    Y = linearmodels(idl).ydata{1,1};
    
    Y_t = Y;
    X_t = X(:,1:4*idm+1);
    
    % constrained linear regression matlab function
    options=optimset;
    options.Algorithm='interior-point'; % to get rid of some warnings
    options.Display='off';
    Aineq = diag([0;repmat([0;-1;1;-1], [idm,1])]);
    bineq = zeros(4*idm+1,1);
    coeff = lsqlin(X_t,Y_t,Aineq,bineq,[],[],[],[],[],options);
    linearmodels(idl).coeff = coeff;
    
    linearmodels(idl).leaves = [];
    linearmodels(idl).xdata = [];
    linearmodels(idl).ydata = [];
        
    progress = idl/numleafs;
    waitbar(progress, h, sprintf('Leaf %d of %d', idl, numleafs));
    
end

close(h);