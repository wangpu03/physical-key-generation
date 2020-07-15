clc;
clear;
addpath(genpath('E:\Github\physical-key-generation\MIToolbox-master\matlab'));
%calculate the entropy of norm h(x)=1/2*log2(2*pi*exp(1)*sigma^2)

%% calculate the entropy of the norm variables with definition
sigma1 = 3;
h_norm1= 1/2*log2(2*pi*exp(1)*sigma1^2);

%calculate the entropy of the norm variables with probability desity function (PDF)
interval = 0.001;
x = -10*sigma1:interval:10*sigma1;
y = normpdf(x, 0, sigma1);
h_inter = sum(y.*log2(1./y))*interval;  %利用积分计算熵

%calculate the entropy of variable with statistics
r1 = normrnd(0,sigma1,[1 10000000])';
% figure(1);
% edges = [-10*sigma:1:10*sigma];
% histogram(r1,edges);
pre = 100;
h_r1 = h(r1*pre)-log2(pre);

%% aother Gaussian variable with variabce 3
sigma2 = 3;
h_norm2= 1/2*log2(2*pi*exp(1)*sigma2^2); % 标准计算熵的公式
r2 = normrnd(0,sigma2,[1 10000000])';
h_r2 = h(r2*pre)-log2(pre);  %利用工具函数

%% joint entropy of two Gaussian variables
r = [r1,r2];
h_joint = h(r);


