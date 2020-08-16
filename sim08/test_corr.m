clc;
clear;
h_ab = randn(1,5000) + 1i*randn(1,5000);
r2_d = 0.9;
h_ae = r2_d*h_ab + sqrt(1-r2_d^2)*(randn(1,5000) + 1i*randn(1,5000));
cor = corrcoef(h_ab,h_ae)
cor_abs = corrcoef(abs(h_ab),abs(h_ae))

exp_PDP = exp(-(0:10-1));
exp_PDP = exp_PDP/norm(exp_PDP);

h_ab_ofdm=(randn(5000,1)+1i*randn(5000,1))*sqrt(exp_PDP/2);

h_ae_ofdm = zeros(5000,10);
for i = 1:length(exp_PDP)
    h_ae_ofdm(:,i) = r2_d*h_ab_ofdm(:,i) + sqrt(exp_PDP(i)/2)*sqrt(1-r2_d^2)*(randn(5000,1) + 1i*randn(5000,1));
end

cor_ofdm = corrcoef(h_ab_ofdm,h_ae_ofdm)


% h_ae_ofdm2 = zeros(5000,10);
% h_ae_ofdm2 = r2_d*h_ab_ofdm .+ sqrt(1-r2_d^2)*randn(1,1) + 1i*randn(1,1)

