%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

%Pt = 110;           %transmission power of the RF source
rho = 0.2;          %power split ratio for information receiver, range is [0,1]

%fig.1 
Pt = 10:20:190;
c_real = zeros(1,10);
c_theo = zeros(1,10);
c_sk = zeros(1,10);
energy = zeros(1,10);
for i = 1:10
    [c_real(i),c_theo(i),c_sk(i), energy(i)] = key_gen_func(rho, Pt(i));
end

yyaxis left
plot(Pt,c_real);
hold on;
plot(Pt,c_theo);
%plot(Pt,c_sk);
yyaxis right
plot(Pt,energy);
hold off

%[c_real,c_theo,c_sk, energy] = key_gen_func(rho, Pt)