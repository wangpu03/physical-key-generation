clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));
%entropy of norm h(x)=1/2*log2(2*pi*exp(1)*sigma^2)

%calculate the entropy of the norm variables with definition
sigma1 = 3;
%h_norm1= 1/2*log2(2*pi*exp(1)*sigma1^2);

%calculate the entropy of the norm variables with probability desity function (PDF)
%interval = 0.001;
%x = -10*sigma1:interval:10*sigma1;
%y = normpdf(x, 0, sigma1);
%h_inter = sum(y.*log2(1./y))*interval;

%calculate the entropy of variable with statistics
r1 = normrnd(0,sigma1,[1 100000000])';
% figure(1);
% edges = [-10*sigma:1:10*sigma];
% histogram(r1,edges);
%pre = 100;
%h_r1 = h(r1*pre)-log2(pre);

sigma2 = 3;
%h_norm2= 1/2*log2(2*pi*exp(1)*sigma2^2);
r2 = normrnd(0,sigma2,[1 100000000])';
%h_r2 = h(r2*pre)-log2(pre);

r = [r1,r2];
h_joint = h(r);



% r_pro = r1.*r2;
% h_r_pro = h(r_pro*pre)-log2(pre);
% cov(r1,r2)
% 
% var(r_pro)
% 
% 
% %==============
% clc;
% clear all;
% f = [119 123 168 119;123 119 168 168;
%      119 119 107 119;107 107 119 119];%?f??????
% p = hist(f(:),8);%???????8?????????????hist(f(:),256)???????256????
% p = p/sum(p);
% i = find(p);
% h = -sum(p(i).*log2(p(i)))%?????

