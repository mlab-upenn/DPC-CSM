function linearmodels = train_linearmodel_in_leaves(model, xTrain_c_all, yTrain_all, idm)

xTrain_d = model.X;
NTrees = model.NTrees;

h0 = waitbar(0, 'Synthesis Forest');

for idt = 1:NTrees
    
    leaf_index = find((model.Trees{idt}.Children(:,1)==0)&(model.Trees{idt}.Children(:,2)==0));
    numleafs = length(leaf_index);
    [~,node] = predict(model.Trees{idt}, xTrain_d(~model.OOBIndices(:,idt),:));
    xTrain_c = xTrain_c_all(~model.OOBIndices(:,idt),:);
    yTrain = yTrain_all(~model.OOBIndices(:,idt),:);
    h = waitbar(0, 'Sythesis Tree');
    
    for idl=1:numleafs
        
        % find indices of samples which end up in this leaf
        leafmodels(idl).leaves = {find(node==leaf_index(idl))}; %#ok<*AGROW>
        
        % mean prediction at this leaf
        leafmodels(idl).mean = model.Trees{idt}.NodeMean(leaf_index(idl));
        
        % the control variables sample values which contribute to this leaf (support)
        leafmodels(idl).xdata = {xTrain_c(leafmodels(idl).leaves{1,1},:)};
        
        % the response variable value which contribute to this leaf
        leafmodels(idl).ydata = {yTrain(leafmodels(idl).leaves{1,1})};
        
        % train a linear model
        X = [ones(size(leafmodels(idl).ydata{1},1),1), leafmodels(idl).xdata{1}];
        Y = leafmodels(idl).ydata{1,1};
        
        Y_t = Y;
        X_t = X(:,1:4*idm+1);
        if size(X_t,1)<size(X_t,2)
            % keyboard;
            leafmodels(idl).coeff = nan(4*idm+1,1);
        else
            %         while size(X_t,1)<size(X_t,2)
            %             ncur = size(X_t,1);
            %             Y_t(ncur+1) = mean(Y_t);
            %             X_t(ncur+1,:) = [1, mean(X_t(:,2:end),1)];
            %         end
            
            % constrained linear regression matlab function
            options=optimset;
            options.Algorithm='interior-point'; % to get rid of some warnings
            options.Display='off';
            Aineq = diag([0;repmat([0;-1;1;-1], [idm,1])]);
            bineq = zeros(4*idm+1,1);
            coeff = lsqlin(X_t,Y_t,Aineq,bineq,[],[],[],[],[],options);
            leafmodels(idl).coeff = coeff;
        end
        
        leafmodels(idl).leaves = [];
        leafmodels(idl).xdata = [];
        leafmodels(idl).ydata = [];
        
        progress = idl/numleafs;
        waitbar(progress, h, sprintf('Leaf %d of %d', idl, numleafs));
        
    end
    
    linearmodels{idt} = leafmodels;
    close(h);
    
    progress = idt/NTrees;
    waitbar(progress, h0, sprintf('Tree %d of %d', idt, NTrees));
    
end

close(h0);

