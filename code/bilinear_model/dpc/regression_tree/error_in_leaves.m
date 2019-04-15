function error_in_leaves(leafmodels, idl, idp)

% extract X,Y from the leaf
X = leafmodels(idl).xdata{1,1};
Y = leafmodels(idl).ydata{1,1};

Y_t = Y;
X_t = X(:,1:idp);
yTrue = zeros(1,size(Y_t,1));
yPred = zeros(1,size(Y_t,1));
for kk = 1:size(Y_t,1)
    yTrue(kk) = Y_t(kk);
    coeff = leafmodels(idl).coeff;
    yPred(kk) = [1, X_t(kk,:)]*coeff;
end

% true-predicted plot after linear regression
figure; hold on; grid on;
h1 = plot(yTrue, 'LineWidth', 1.5);
h2 = plot(yPred, 'LineWidth', 1.5);
legend([h1, h2], 'true', 'predicted');




