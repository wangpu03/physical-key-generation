%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

Pt = 20;           %transmission power of the RF source
%rho = 0.2;          %power split ratio for information receiver, range is [0,1]

%fig.1 
rho = 0:0.04:1;
c_real = zeros(1,length(rho));
c_theo = zeros(1,length(rho));
c_sk = zeros(1,length(rho));
energy = zeros(1,length(rho));
for i = 1:length(rho)
    [c_real(i),c_theo(i),c_sk(i), energy(i)] = key_gen_func(rho(i), Pt);
end

figure(1);
yyaxis left
plot(rho,c_real,'-s');
hold on;
plot(rho,c_theo,'-<');
%plot(rho,c_sk,'->');
yyaxis right
plot(rho,energy,'-o');
hold off
grid on;

lambda1 = 1;
G1 = (1-lambda1)*(c_real)/max(c_real)+lambda1*(energy/max(energy));

lambda2 = 0.6;
G2 = (1-lambda2)*(c_real)/max(c_real)+lambda2*(energy/max(energy));

lambda3 = 0.5;
G3 = (1-lambda3)*(c_real)/max(c_real)+lambda3*(energy/max(energy));

lambda4 = 0.4;
G4 = (1-lambda4)*(c_real)/max(c_real)+lambda4*(energy/max(energy));

lambda5 = 0;
G5 = (1-lambda5)*(c_real)/max(c_real)+lambda5*(energy/max(energy));
figure(2);
plot(rho,G1,'-o');
hold on;
plot(rho,G2,'-+');
hold on;
plot(rho,G3,'-*');
plot(rho,G4,'-d');
grid on;

plot(rho,G5,'-s');
hold off
