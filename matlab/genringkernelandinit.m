% generate data points (rings)
rng(200);
randdata = rand(10000,2);
theta = pi*rand(300,1);
distances = rand(300,1)*4+7;
X = [distances.*cos(theta) distances.*sin(theta)];
theta = pi*rand(300,1);
distances = rand(300,1)*4;
X = [X; distances.*cos(theta) distances.*sin(theta)];
N = size(X,1);

% random initial clustering
clustering = randi([1,2],N,1);
% visualize inital clustering
figure;plot(X(clustering==2,1),X(clustering==2,2),['r','s'], 'MarkerFaceColor','r','MarkerSize',5); hold on; 
plot(X(clustering==1,1),X(clustering==1,2),['b','s'], 'MarkerFaceColor','b','MarkerSize',5); axis equal
title('initial clustering');

% Gaussian affinity matrix
A = zeros(N,N,'double');
sigma = 3;
for n1=1:N
   xn1 = X(n1,:);
   diff = repmat(xn1,N,1)-X;
   A(n1,:) = exp(-sum(diff.^2,2)/2/sigma/sigma);
end