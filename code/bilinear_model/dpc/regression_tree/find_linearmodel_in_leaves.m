function coeff = find_linearmodel_in_leaves(model, leafmodels, features_d)
    
outputs = predict(model, features_d);

% find the corresponding leaf index
for idl = 1:length(leafmodels)
    if abs(outputs - leafmodels(idl).mean)<1e-8
        break;
    end
end

% litmus test
% [testypred,testleaf] = predict(model,features_d);
% [~, testnode] = predict(model,model.X);
% if mean(model.Y(find(testnode==testleaf)))-testypred>1e-10
%     keyboard;
% end
% if find(testnode==testleaf)~=leafmodels(idl).leaves{1,1}
%     keyboard;
% end

% extract the linear model coefficients from that leaf
coeff = leafmodels(idl).coeff;