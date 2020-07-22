clc;
clear;

addpath('E:\Github\physical-key-generation\MIToolbox-master\matlab');

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

snr = 10:3:40;
v1 = zeros(50000, length(snr));
v2 = zeros(50000, length(snr));
%%
% nfft = 128, K =1
n_fft = 128;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);


for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1128 = v1;
v2_1128 = v2;

v1_bit_1128 = v1_bit;
v2_bit_1128 = v2_bit;
mi_bit_1128 = mi_bit;

%%
% nfft = 128, K =3
n_fft = 128;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 3;
c_flag = 1;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_3128 = v1;
v2_3128 = v2;

v1_bit_3128 = v1_bit;
v2_bit_3128 = v2_bit;
mi_bit_3128 = mi_bit;

%%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1256 = v1;
v2_1256 = v2;

v1_bit_1256 = v1_bit;
v2_bit_1256 = v2_bit;
mi_bit_1256 = mi_bit;

%%
% nfft = 256, K =3
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=3;
c_flag = 1;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_3256 = v1;
v2_3256 = v2;

v1_bit_3256 = v1_bit;
v2_bit_3256 = v2_bit;
mi_bit_3256 = mi_bit;


%%
% nfft = 128, K = 1， c_flag = 0
n_fft = 128;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 1;
c_flag = 0;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

% for k = 1:50000
%     [v11(k),v22(k)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, 50);
% end 

for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1128_guard = v1;
v2_1128_guard = v2;

v1_bit_1128_guard = v1_bit;
v2_bit_1128_guard = v2_bit;
mi_bit_1128_guard = mi_bit;


%%
% nfft = 256, K = 1， c_flag = 0
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 1;
c_flag = 0;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

% for k = 1:50000
%     [v11(k),v22(k)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, 50);
% end 

for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1256_guard = v1;
v2_1256_guard = v2;

v1_bit_1256_guard = v1_bit;
v2_bit_1256_guard = v2_bit;
mi_bit_1256_guard = mi_bit;


%%
% nfft = 128, K = 1， c_flag = 0， noiseless
n_fft = 128;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 0;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);

snr = 2000:300:5000;
for i = 1:length(snr)
for j = 1:50000
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), 50000,1);
mean_v2 = repmat(mean(v2), 50000,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1128_guard_n = v1;
v2_1128_guard_n = v2;

v1_bit_1128_guard_n = v1_bit;
v2_bit_1128_guard_n = v2_bit;
mi_bit_1128_guard_n = mi_bit;


% %%
% % nfft = 256, K = 3， c_flag = 0， noiseless
% n_fft = 256;
% n_cp = n_fft/4;    % size of cyclic prefix extension
% n_ofdm = n_fft + n_cp;
% K=1;
% c_flag = 0;
% % generate the data 
% rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
% save data_input.txt -ascii rand_ints_gen
% rand_ints = load("data_input.txt");
% data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);
% 
% for i = 1:length(snr)
% for j = 1:50000
%     [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
% end
% end
% 
% mean_v1 = repmat(mean(v1), 50000,1);
% mean_v2 = repmat(mean(v2), 50000,1);
% 
% v1_bit = double(v1>= mean_v1);
% v2_bit = double(v2>= mean_v2);
% 
% for i = 1:length(snr)
%     mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
% end
% 
% v1_3256_guard_n = v1;
% v2_3256_guard_n = v2;
% 
% v1_bit_3256_guard_n = v1_bit;
% v2_bit_3256_guard_n = v2_bit;
% mi_bit_3256_guard_n = mi_bit;

%%
snr = 10:3:40;
plot(snr, mi_bit_1128,'k-o','LineWidth',1);
hold on;
plot(snr, mi_bit_3128,'k--x','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1256,'r-d','LineWidth',1.5);
hold on;
plot(snr, mi_bit_3256,'r--s','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1128_guard,'b-^','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1256_guard,'b--v','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1128_guard_n(:,1:11),'y--*','LineWidth',1.5);
grid on;
% plot(snr, mi_bit_3256_guard_n(:,1:11),'y--*','LineWidth',1);
%axis([20 50 0.4 1]);
legend('K=1, N_{FFT} = 128','K=3, N_{FFT} = 128','K=1, N_{FFT} = 256','K=3, N_{FFT} = 256','K=1, N_{Guard} = 128','K=1, N_{Guard} = 256','K=1, N_{Guard}= 128,Noiseless');
ylabel('Mutual Information per bit');
xlabel('SNR (dB)');

%% 计算BDR
bdr_1128 = 1 - sum(v1_bit_1128 == v2_bit_1128)/50000;
bdr_3128 = 1 - sum(v1_bit_3128 == v2_bit_3128)/50000;
bdr_1256 = 1 - sum(v1_bit_1256 == v2_bit_1256)/50000;
bdr_3256 = 1 - sum(v1_bit_3256 == v2_bit_3256)/50000;
bdr_1128_guard = 1 - sum(v1_bit_1128_guard == v2_bit_1128_guard)/50000;
bdr_1256_guard = 1 - sum(v1_bit_1256_guard == v2_bit_1256_guard)/50000;
bdr_1128_guard_n = 1 - sum(v1_bit_1128_guard_n == v2_bit_1128_guard_n)/50000;

plot(snr, bdr_1128,'k-o','LineWidth',1);
hold on;
plot(snr, bdr_3128,'k--x','LineWidth',1);
hold on;
plot(snr, bdr_1256,'r-d','LineWidth',1);
hold on;
plot(snr, bdr_3256,'r--s','LineWidth',1);
hold on;
plot(snr, bdr_1128_guard,'b-^','LineWidth',1);
hold on;
plot(snr, bdr_1256_guard,'b--v','LineWidth',1);
hold on;
plot(snr, bdr_1128_guard_n,'y--*','LineWidth',1);
grid on;
%axis([20 50 0.6 1]);
legend('K=1, N_{FFT} = 128','K=3, N_{FFT} = 128','K=1, N_{FFT} = 256','K=3, N_{FFT} = 256','K=1, N_{Guard} = 128','K=1, N_{Guard} = 256','K=1, N_{Guard}= 128,Noiseless');
ylabel('Bits Disagreement Ratio');
xlabel('SNR (dB)');
