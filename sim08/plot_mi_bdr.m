clc;
clear;
datestr(now)
addpath('E:\Github\physical-key-generation\MIToolbox-master\matlab');

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

% 该数据可以生成一次之后可以永久使用，保证每次调制的数据一致
% rand_ints_gen = randi(2,10000,1)-1;
% save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");

num_mc = 50000;

snr = 0:3:30;
v1 = zeros(num_mc, length(snr));
v2 = zeros(num_mc, length(snr));
%%
% nfft = 64, K =1
n_fft = 64;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_164 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_164, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm, n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_164 = v1;
v2_164 = v2;

v1_bit_164 = v1_bit;
v2_bit_164 = v2_bit;
mi_bit_164 = mi_bit;

%%
% nfft = 64, K =3
n_fft = 64;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 3;
c_flag = 1;
% generate the data 
rand_ints_364 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_364, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_364 = v1;
v2_364 = v2;

v1_bit_364 = v1_bit;
v2_bit_364 = v2_bit;
mi_bit_364 = mi_bit;
plot(snr, mi_bit_164,'k-o','LineWidth',1.5);
hold on;
plot(snr, mi_bit_364,'k--x','LineWidth',1.5);

%%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

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
rand_ints_3256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_3256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

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


%% ===========================================
% nfft = 64, K =1, low multipath spread
n_fft = 64;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_164 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_164, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr2(data_ofdm, n_ofdm, n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_164_2 = v1;
v2_164_2 = v2;
v1_bit_164_2 = v1_bit;
v2_bit_164_2 = v2_bit;
mi_bit_164_2 = mi_bit;

%%
% nfft = 64, K =3,low multipath spread
n_fft = 64;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 3;
c_flag = 1;
% generate the data 
rand_ints_364 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_364, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr2(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_364_2 = v1;
v2_364_2 = v2;

v1_bit_364_2 = v1_bit;
v2_bit_364_2 = v2_bit;
mi_bit_364_2 = mi_bit;

%%
% nfft = 256, K =1,low multipath spread
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=1;
c_flag = 1;
% generate the data 
rand_ints_1256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_1256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr2(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_1256_2 = v1;
v2_1256_2 = v2;
v1_bit_1256_2 = v1_bit;
v2_bit_1256_2 = v2_bit;
mi_bit_1256_2 = mi_bit;

%%
% nfft = 256, K =3, low multipath spread
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K=3;
c_flag = 1;
% generate the data 
rand_ints_3256 = rand_ints(1:n_fft*mod_order*K,:);
data_ofdm = ofdm_module(rand_ints_3256, mod_method, n_fft, n_cp, c_flag);

for i = 1:length(snr)
for j = 1:num_mc
    [v1(j,i),v2(j,i)] = mi_snr2(data_ofdm, n_ofdm,n_cp, K, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_mc,1);
mean_v2 = repmat(mean(v2), num_mc,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);

for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

v1_3256_2 = v1;
v2_3256_2 = v2;
v1_bit_3256_2 = v1_bit;
v2_bit_3256_2 = v2_bit;
mi_bit_3256_2 = mi_bit;


%%
snr = 0:3:30;
plot(snr, mi_bit_164,'k-o','LineWidth',1.5);
hold on;
plot(snr, mi_bit_364,'r-x','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1256,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, mi_bit_3256,'-s','LineWidth',1.5,'Color', '#0072BD');
hold on;
plot(snr, mi_bit_164_2,'k--o','LineWidth',1.5);
hold on;
plot(snr, mi_bit_364_2,'r--x','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1256_2,'--d','LineWidth',1.5,'Color', '#77AC30'); %+0.03./(snr./3+1)
hold on;
plot(snr, mi_bit_3256_2,'--s','LineWidth',1.5,'Color', '#0072BD'); %+0.02./(snr./3+1)
grid on;
axis([0 30 0.5 0.95]);
legend('K = 1, N = 64, L = 10','K = 3, N = 64, L = 10','K = 1, N = 256, L = 10','K = 3, N = 256, L = 10','K = 1, N = 64, L = 5','K = 3, N = 64, L = 5','K = 1, N = 256, L = 5','K = 3, N = 256, L = 5');
ylabel('Mutual Information per bit');
xlabel('SNR (dB)');

%% 计算BDR
bdr_164 = 1 - sum(v1_bit_164 == v2_bit_164)/num_mc;
bdr_364 = 1 - sum(v1_bit_364 == v2_bit_364)/num_mc;
bdr_1256 = 1 - sum(v1_bit_1256 == v2_bit_1256)/num_mc;
bdr_3256 = 1 - sum(v1_bit_3256 == v2_bit_3256)/num_mc;
bdr_164_2 = 1 - sum(v1_bit_164_2 == v2_bit_164_2)/num_mc;
bdr_364_2 = 1 - sum(v1_bit_364_2 == v2_bit_364_2)/num_mc;
bdr_1256_2 = 1 - sum(v1_bit_1256_2 == v2_bit_1256_2)/num_mc;
bdr_3256_2 = 1 - sum(v1_bit_3256_2 == v2_bit_3256_2)/num_mc;

%%
snr = 0:3:30;
plot(snr, bdr_164,'k-o','LineWidth',1.5);
hold on;
plot(snr, bdr_364,'r-x','LineWidth',1.5);
hold on;
plot(snr, bdr_1256,'-d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, bdr_3256,'-s','LineWidth',1.5,'Color', '#0072BD');
hold on;
plot(snr, bdr_164_2,'k--o','LineWidth',1.5);
hold on;
plot(snr, bdr_364_2,'r--x','LineWidth',1.5);
hold on;
plot(snr, bdr_1256_2,'--d','LineWidth',1.5,'Color', '#77AC30');
hold on;
plot(snr, bdr_3256_2,'--s','LineWidth',1.5,'Color', '#0072BD');
grid on;
axis([0 30 0.0 0.08]);
legend('K = 1, N = 64, L = 10','K = 3, N = 64, L = 10','K = 1, N = 256, L = 10','K = 3, N = 256, L = 10','K = 1, N = 64, L = 5','K = 3, N = 64, L = 5','K = 1, N = 256, L = 5','K = 3, N = 256, L = 5');
ylabel('Bit Disagreement Ratio');
xlabel('SNR (dB)');
datestr(now)