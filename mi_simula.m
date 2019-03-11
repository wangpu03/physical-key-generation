%% the simulation about key generation between two backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

sigma_f = 4;
f = normrnd(0,sigma_f,[1 100000000])';

sigma_g = 3;
g = normrnd(0,sigma_f,[1 100000000])';

sigma_za = 3;
za = normrnd(0,sigma_za,[1 100000000])';

sigma_zb = 3;
zb = normrnd(0,sigma_zb,[1 100000000])';

X = f.*g + za;
Y = f.*g + zb;

pre = 100;
h_x = h(X*pre)-log2(pre);
h_theo_x = 1/2*log2(2*pi*exp(1)*(var(f.*g)+var(za)))
h_y = h(Y*pre)-log2(pre)

mi_theo= 1/2*log2(1+(var(f.*g)/(var(za)+var(zb)+(var(za)*var(zb)/var(f.*g)))))

mi_xy = mi(X,Y);

