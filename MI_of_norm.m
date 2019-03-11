%% the simulation about the basic model of physical layer key generation scheme
% Two independent complex Gaussian random variables Z1~N(0,N1),Zb~N(0,N2),
% serving as noises corrupting the channel observation by Alice and Bob,
% respectively.
% A complex Guassian random variables H~N(0,p) as the joint channel observation 
clc;
clear;
addpath(genpath('C:\Code file\Matlab\keygen\MIToolbox-master\matlab'));
%entropy of norm h(x)=1/2*log2(2*pi*exp(1)*sigma^2)

%calculate the entropy of the norm variables with definition
% signal n1
sigma1 = 3;
h_r1_theo= 1/2*log2(2*pi*exp(1)*sigma1^2)
%calculate the entropy of variable with statistics
r1 = normrnd(0,sigma1,[1 100000000])';
pre = 100;
h_r1 = h(r1*pre)-log2(pre)
%h_r11 = h(r1)

% signal n2
sigma2 = 3;
h_r2_theo= 1/2*log2(2*pi*exp(1)*sigma2^2)
r2 = normrnd(0,sigma2,[1 100000000])';
h_r2 = h(r2*pre)-log2(pre)

% joint entropy of signal n1 and signal n2
r = [r1,r2];
hj_r1r2 = h(r)

%%
% signal H
sigma_h = 30;
r_h = normrnd(0,sigma_h,[1 100000000])';

%received total signal x and y
x = r_h+r1;
y = r_h+r2;
h_x = h(x*pre)-log2(pre); %the entropy of x
h_y = h(y*pre)-log2(pre); %the entropy of y
h_x_theo= 1/2*log2(2*pi*exp(1)*(sigma_h^2+sigma1^2)) %the theoretical entropy of x

hj_xy =h([x,y]);% h([x*pre,y*pre])-2*log2(pre) %the joint entropy of x and y
%the theoretical joint entropy of x and y
hj_xy_theo = 1/2*log2((2*pi*exp(1))^2*((sigma_h^2+sigma1^2)*(sigma_h^2+sigma2^2)-sigma_h^4));

%% the mutual informatino between x and y
%  calculate from the I(X;Y) = h(x)+h(Y)-H(X,Y)
%  calculate derictly with the function
%  calculate with the theorectical formulation 
mi_cal_xy = h_x+h_y-hj_xy
mi_xy = mi(x,y) %calculate derictly from the function
mi_theo = 1/2*log2(1+(sigma_h^2/(sigma1^2+sigma2^2+(sigma1^2*sigma2^2/sigma_h^2))))


