% Identify numeric data for PCA
X = data{:,3:end};

%% Perform PCA coordinate transformation
[pcs,scrs,vexp,~,pexp] = pca(X);

%% Visualize variance explained by component
figure
pareto(pexp)
xlabel('Principal Component')
ylabel('Percent Variance Explained')
title('PCA Percent Variance Explained')

%% Visualize first 3 principal components with biplot
figure
subplot(2,2,1)
biplot(pcs(:,1:3)) 
xlabel('Component 1')
ylabel('Component 2')
zlabel('Component 3')

subplot(2,2,2)
biplot(pcs(:,[1,3]))
xlabel('Component 1')
ylabel('Component 3')

subplot(2,2,3)
biplot(pcs(:,[2,3])) 
xlabel('Component 2')
ylabel('Component 3')

subplot(2,2,4)
biplot(pcs(:,[1,2])) 
xlabel('Component 1')
ylabel('Component 2')

%% Visualize with variable names on plot
% Customize the biplot visualization
figure
h = biplot(pcs(:,1:3),'VarLabels',features);
title('First 3 Principal Components')

hlines = h(1:42);
hmarkers = h(43:84);
hlabels = h(85:end);

% Keep only variable names with Grav* in the name
idx = contains(features,'Grav');
delete(hlabels(~idx))

% Change marker properties 
set(hmarkers,'Color','r','Marker','*')

%% Visualize pca coefficients of each variable with heatmap
f = figure;
f.Position = f.Position.*[1 1 1 2];
imagesc(abs(pcs(:,1:3)))
colorbar
yticks(1:42)
yticklabels(features)
xticks(1:3)
xlabel('Principal Component')
ylabel('Variable')
title('First Three Principal Components')

%% Visualize pca scores with parallel coordinates plot
figure
parallelcoords(scrs,'Group',data.Activity,'Quantile',0.25)
xticklabels([])
title('PCA Scores Per Variable')

% You can also visualize pca scores directly from original data with parallel coordinates plot
% parallelcoords(X,'Group',data.Activity,'Standardize','PCA','Quantile',0.25)
