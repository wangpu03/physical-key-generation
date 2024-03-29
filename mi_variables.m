%% the simulation about key generation between two backscatter devices
% calculate the mutual information between two variables X, Y
% f channel information (Guassin variables)
% g channel information (Guassin variables)
% z1 and z1 noise signals ((Guassin variables))
% X = f.*g + z1;
% Y = f.*g + z2;
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

sigma_f = 7;
f = normrnd(0,sigma_f,[1 100000000])';

sigma_g = 6;
g = normrnd(0,sigma_g,[1 100000000])';

sigma_z1 = 1;
z1 = normrnd(0,sigma_z1,[1 100000000])';

sigma_z2 = 1;
z2 = normrnd(0,sigma_z2,[1 100000000])';

sigma_z3 = 1;
z2 = normrnd(0,sigma_z2,[1 100000000])';

X = f.*g + z1;
Y = f.*g + z2;

pre = 100;
h_x = h(X*pre)-log2(pre);
h_theo_x = 1/2*log2(2*pi*exp(1)*((sigma_f^2*sigma_g^2)+sigma_z1^2));
h_y = h(Y*pre)-log2(pre);
h_joint_xy = h([X,Y]);

mi_theo= 1/2*log2(1+(var(f.*g)/(var(z1)+var(z2)+(var(z1)*var(z2)/var(f.*g)))));

mi_xy = mi(X,Y);

