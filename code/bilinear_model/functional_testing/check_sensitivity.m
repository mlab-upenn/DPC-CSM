
u = zeros(4,length(time));
x = zeros(12,length(time));
y = zeros(3,length(time));
x(:,1) = [22 22 22 22 22 22 22 22 22 20 15 10]';
y(:,1) = [22 22 22]';

for k = 1:length(time)-1
[x(:,k+1), y(:,k+1)] = simulate_model(model, x(:,k), d(:,k), u(:,k));
end
figure;
plot(1:200,y(1,1:200))