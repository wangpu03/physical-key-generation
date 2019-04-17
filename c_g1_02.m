%% simulation of the key generation between two Backscatter devices
clc;
clear;
addpath(genpath('C:\Code file\Matlab\physical-key-generation\MIToolbox-master\matlab'));

log_pt = [30,33];
Pt = 10.^(log_pt/10);                %transmission power of the RF source
%rho = 0.2;             %power split ratio for information receiver, range is [0,1]

d_f1 = 5.1;
d_f2 = 6;
d_h12 = 3;

rho = 0:0.1:1;
c_real = zeros(2,length(rho));
c_theo = zeros(2,length(rho));
e_real = zeros(2,length(rho));
e_theo = zeros(2,length(rho));

for j = 1:1:2
    for i = 1:length(rho)
        [c_real(j,i),c_theo(j,i),e_real(j,i), e_theo(j,i)] = key_gen_func_02(d_f1,d_f2,d_h12,rho(i), Pt(j));
    end
end

G1 = zeros(2,length(rho));
G3 = zeros(2,length(rho));
G5 = zeros(2,length(rho));
% need to run it gain and plot with different SNR
for j = 1:1:2
    for i = 1:length(rho)
        lambda1 = 0.7;
        G1(j,i) = (1-lambda1)*(c_real(j,i))/max(c_real(2,:))+lambda1*(e_real(j,i)/max(e_real(2,:)));
        lambda3 = 0.4;
        G3(j,i) = (1-lambda3)*(c_real(j,i))/max(c_real(2,:))+lambda3*(e_real(j,i)/max(e_real(2,:)));
        lambda5 = 0.1;
        G5(j,i) = (1-lambda5)*(c_real(j,i))/max(c_real(2,:))+lambda5*(e_real(j,i)/max(e_real(2,:)));
    end
end

figure(1);
plot(rho,G1(1,:),'-o','DisplayName','\lambda = 0.7,SNR = 30','LineWidth',1);
hold on;
plot(rho,G3(1,:),'-*','DisplayName','\lambda = 0.5,SNR = 30','LineWidth',1);
plot(rho,G5(1,:),'-d','DisplayName','\lambda = 0.3,SNR = 30','LineWidth',1);

plot(rho,G1(2,:),'--o','DisplayName','\lambda = 0.7,SNR = 33','LineWidth',1);
plot(rho,G3(2,:),'--*','DisplayName','\lambda = 0.5,SNR = 33','LineWidth',1);
plot(rho,G5(2,:),'--d','DisplayName','\lambda = 0.3,SNR = 33','LineWidth',1);

grid on;
hold off;
xlabel('\rho');
ylabel ('Objective G1');
legend;