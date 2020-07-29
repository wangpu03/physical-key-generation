clc;
clear;

addpath('E:\Github\physical-key-generation\MIToolbox-master\matlab');

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));

%
% nfft = 256, K =1
n_fft = 256;
n_cp = n_fft/4;    % size of cyclic prefix extension
n_ofdm = n_fft + n_cp;
K = 1;
c_flag = 1;
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*K,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");
data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, c_flag);
num_sim = 50000;
snr = 10:3:40;
v1 = zeros(num_sim, length(snr));
v2 = zeros(num_sim, length(snr));
ve1 = zeros(num_sim, length(snr));
ve2 = zeros(num_sim, length(snr));

%%
d = [8 7 3];
taps = [8 7 3];
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = mi_dist(data_ofdm, n_ofdm, n_cp ,K, d, taps, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_873 = v1;
v2_873 = v2;
ve1_873 = ve1;
ve2_873 = ve2;
mi_bit_873 = mi_bit;
cmi_ve2_873 = cmi_ve2;
cmi_ve1_873 = cmi_ve1;

%%
d = [12 7 3];
taps = [12 7 3];
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = mi_dist(data_ofdm, n_ofdm, n_cp ,K, d, taps, snr(i));
end
end
mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_1273 = v1;
v2_1273 = v2;
ve1_1273 = ve1;
ve2_1273 = ve2;
mi_bit_1273 = mi_bit;
cmi_ve2_1273 = cmi_ve2;
cmi_ve1_1273 = cmi_ve1;


%%
d = [12 7 6];
taps = [12 7 6];
for i = 1:length(snr)
for j = 1:num_sim
    [v1(j,i),v2(j,i),ve1(j,i),ve2(j,i)] = mi_dist(data_ofdm, n_ofdm, n_cp ,K, d, taps, snr(i));
end
end

mean_v1 = repmat(mean(v1), num_sim,1);
mean_v2 = repmat(mean(v2), num_sim,1);
mean_ve1 = repmat(mean(ve1), num_sim,1);
mean_ve2 = repmat(mean(ve2), num_sim,1);

v1_bit = double(v1>= mean_v1);
v2_bit = double(v2>= mean_v2);
ve1_bit = double(ve1>= mean_ve1);
ve2_bit = double(ve2>= mean_ve2);

cmi_ve2 = zeros(1,length(snr));
cmi_ve1 = zeros(1,length(snr));
mi_bit = zeros(1,length(snr));
for i = 1:length(snr)
    mi_bit(i) = mi(v1_bit(:,i),v2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve2(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve2_bit(:,i));
end

for i = 1:length(snr)
    cmi_ve1(i) = cmi(v1_bit(:,i),v2_bit(:,i),ve1_bit(:,i));
end

v1_1276 = v1;
v2_1276 = v2;
ve1_1276 = ve1;
ve2_1276 = ve2;
mi_bit_1276 = mi_bit;
cmi_ve2_1276 = cmi_ve2;
cmi_ve1_1276 = cmi_ve1;

plot(snr, mi_bit_873,'r-o','LineWidth',1.5);
hold on;
plot(snr, mi_bit_1273,'k-x','LineWidth',1.5);
plot(snr, mi_bit_1276,'b-*','LineWidth',1.5);
plot(snr, cmi_ve2_873,'r--s','LineWidth',1.5);
plot(snr, cmi_ve2_1273,'k--^','LineWidth',1.5);
plot(snr, cmi_ve2_1276,'b--v','LineWidth',1.5);
plot(snr, mi_bit_873-cmi_ve2_873,'r-.v','LineWidth',1.5);
plot(snr, mi_bit_1273-cmi_ve2_1273,'k-.s','LineWidth',1.5);
plot(snr, mi_bit_1276-cmi_ve2_1276,'b-.d','LineWidth',1.5);
hold off;
grid on;

axis([10 40 0 1.1])
txt = {'Red: d_1 = 8, d_2 = 7, d_{12} = 3 ','Black: d_1 = 12, d_2 = 7, d_{12} = 3','Blue: d_1 = 12, d_2 = 7, d_{12} = 6'};
annotation('textbox',[0.15,0.74,0.35,0.165],'LineStyle','-','LineWidth',1,'String',txt)

elps1 = annotation('ellipse',[0.800 0.715 0.0500 0.0500]);
elps1.LineWidth = 1.5;
elps1.Color = '#D95319';

x = [0.770,0.820];
y = [0.57,0.715];
a1 = annotation('textarrow',x, y,'String','Secret key rate');
a1.Color = 'red';
a1.FontSize = 12;
a1.LineWidth = 1.5;

elps2 = annotation('ellipse',[0.6500 0.75 0.0500 0.0500]);
elps2.LineWidth = 1.5;
elps2.Color = '#D95319';

x2 = [0.65,0.6700];
y2 = [0.65,0.75];
a2 = annotation('textarrow',x2, y2,'String','Mutual information');
a2.Color = 'red';
a2.FontSize = 12;
a2.LineWidth = 1.5;

elps3 = annotation('ellipse',[0.800 0.117 0.0500 0.0500]);
elps3.LineWidth = 1.5;
elps3.Color = '#D95319';
x3 = [0.77,0.830];
y3 = [0.25,0.17];
a3 = annotation('textarrow',x3, y3,'String','Leak information');
a3.Color = 'red';
a3.FontSize = 12;
a3.LineWidth = 1.5;
%ylabel('Information (bit)');
xlabel('SNR (dB)');

