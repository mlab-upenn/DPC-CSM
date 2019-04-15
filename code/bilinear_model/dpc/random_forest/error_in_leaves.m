function error_in_leaves(leafmodels, idl, idp)

% extract X,Y from the leaf
X = leafmodels(idl).xdata{1,1};
X = [ones(size(leafmodels(idl).ydata{1},2),1), X'];
Y = leafmodels(idl).ydata{1,1};

for kk = idp %:length(leafmodels(1).mdl)
    yTrue = Y(kk,:)';
    coeff = leafmodels(idl).mdl(kk).coeff;
    yPred = X*coeff;
end

% true-predicted plot after linear regression
figure; hold on; grid on;
h1 = plot(yTrue, 'LineWidth', 1.5);
h2 = plot(yPred, 'LineWidth', 1.5);
legend([h1, h2], 'true', 'predicted');




