clear all; close all; clc;
%%
image_number = 4;
%%
X = imread(['I',num2str(image_number),'_UAV.jpg']); 
Y = imread(['I',num2str(image_number),'_MAN.jpg']); 
%%
x = reshape(im2double(X),[],3);
y = double(~reshape(Y > 220,[],1));
%%
fname = ['Results_k3_I',num2str(image_number)]; % folder name is 'Results_k{number_of_classes}_I{image_number}'
storage_path = fullfile(pwd, fname);
mkdir(fname);
%%
figure(1);
image(X); axis equal;
print(fullfile(storage_path, 'UAV.png'),'-dpng')
print(fullfile(storage_path, 'UAV.eps'),'-depsc')

figure(2);
imagesc((Y > 200)); axis equal;
colormap('gray');
print(fullfile(storage_path, 'MAN.png'),'-dpng')
print(fullfile(storage_path, 'MAN.eps'),'-depsc')
%% GMM
k = 3; % number of classes
T = 5; % number of iterations
tic
[GMM_Mean, GMM_vcov, MODE_conf, y_MODE, mu_hat, pvec] = mygmm(x,y,k,T);
toc
%% MAP Processing
y_MAP3 = getz(pvec);
% G is the panicle class
[~,G] = max(mean([GMM_Mean(1:3) GMM_Mean(4:6) GMM_Mean(7:9)])); %% change according to k
%%
[pos_0,~] = find(y_MAP3 ~= (G-1));
[pos_1,~] = find(y_MAP3 == (G-1));
y_MAP(pos_0) = 0;
y_MAP(pos_1) = 1;
%% MAP Results
final = confusionmat(y,y_MAP);
MAP_conf = 100*final./sum(final,2);
PaddyRatio_Base = 100*sum(y)/length(y);
PaddyRatio_MAP = 100*sum(y_MAP)/length(y);
Pixel_Base = [length(y)-sum(y); sum(y)];
Pixel_MAP = [length(y_MAP)-sum(y_MAP); sum(y_MAP)];
%% MAP Image
figure(3); 
y_MAP_image = reshape(y_MAP,size(X,1),size(X,2));
colormap('gray');
imagesc(mat2gray(~y_MAP_image),[0 1]); axis equal;

print(fullfile(storage_path, 'MAP.png'),'-dpng')
print(fullfile(storage_path, 'MAP.eps'),'-depsc')
%% P-vector : 2 Class
p_1 = pvec(:,G); % 1
p_0 = 1-p_1; % 0
%% ROC Curve
figure(5); 
[X2,Y2,T2,AUC2] = perfcurve(y, p_1, 1);
plot(X2,Y2,'LineWidth',2); hold on;
plot(0:1,0:1,'k','LineWidth',2);
xlabel('False Positive'); ylabel('True Positive');
title('Class: Panicle Pixel (1)');

% print(fullfile(storage_path, 'ROC.png'),'-dpng')
% print(fullfile(storage_path, 'ROC.eps'),'-depsc')
%% p1th Result
y_p1_th = p_1 > 0.999;
p1th = confusionmat(y,double(y_p1_th));
p1th_conf = 100*p1th./sum(p1th,2);
PaddyRatio_p1th = 100*sum(y_p1_th)/length(y);
%% P1th Image
figure(4); 
y_p1th_image = reshape(y_p1_th,size(X,1),size(X,2)); 
colormap('gray');
imagesc(~mat2gray(y_p1th_image),[0 1]); axis equal;

print(fullfile(storage_path, 'p1th.png'),'-dpng')
print(fullfile(storage_path, 'p1th.eps'),'-depsc')
%%
PaddyRatio_Base
PaddyRatio_MAP 
PaddyRatio_p1th
%%
MAP_conf
p1th_conf 
%%
save(fullfile(storage_path,['WS_Image_',num2str(image_number),'.mat']))
%%




