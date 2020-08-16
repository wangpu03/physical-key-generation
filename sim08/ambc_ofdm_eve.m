% consists of all precedures of the OFDM 
% 
clear;
clc;

% simulation parameters

n_fft = 128;        % IFFT/FFT size
n_cp = n_fft/4;    % size of cyclic prefix extension
n_frame = 1;        % the acount of the OFDM frame
n_ofdm = n_fft + n_cp;

d = [8, 7, 3];
taps = [8, 7, 3];
n_L = max(taps(1),taps(2))+taps(3)-1;

h1 = ray_model(d(1),taps(1));
h2 = ray_model(d(2),taps(2));

h_12 = ray_model(d(3),taps(3));
h_21 = h_12;

cor_h1_he = 1;
cor_h12_he2 = 1;

% assume Eve locate closer to device 1
cor = 0;
he = eave_ray_model(h1,cor_h1_he,d(1));
h_e2 = eave_ray_model(h_12,cor_h12_he2,d(3));

d_e1 = 1;
tap_e1 = 1;
h_e1 = ray_model(d_e1,tap_e1);

% modulation methods: BPSK, QPSK,16QAM, 32QAM,64QAM
mod_method = 'QPSK';
% calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK','8PSK','16QAM', '32QAM','64QAM'};
mod_order = find(ismember(mod_methods, mod_method));
% generate the data 
rand_ints_gen = randi(2,n_fft*mod_order*n_frame,1)-1;
save data_input.txt -ascii rand_ints_gen
rand_ints = load("data_input.txt");

data_ofdm = ofdm_module(rand_ints, mod_method, n_fft, n_cp, 0);

pt = 1;
snr = 40;
alpha = 0.3 + 1i*0.4;
pw_noise = pt/10^(snr/10);

% all transmissions in ambient backscatter communication
y1_d = ofdm_trans(data_ofdm, h1, pw_noise);
back1 = ofdm_back(y1_d, alpha, n_ofdm, n_frame);

y2_d = ofdm_trans(data_ofdm, h2, pw_noise);
back2 = ofdm_back(y2_d, alpha, n_ofdm, n_frame);

y1_b = ofdm_trans(back2, h_12, pw_noise);
y2_b = ofdm_trans(back1, h_12, pw_noise);

[y1,y1_cp] = receiver_design(y1_d, y1_b, n_ofdm, n_cp, n_frame, n_L);
[y2,y2_cp] = receiver_design(y2_d, y2_b, n_ofdm, n_cp, n_frame, n_L);

v_y1_cp = mean(abs(y1_cp));
v_y2_cp = mean(abs(y2_cp));

ye_d = ofdm_trans(data_ofdm, he, pw_noise);
y_b_e2 = ofdm_trans(back2, h_e2, pw_noise);
y_b_e1 = ofdm_trans(back1, h_e1, pw_noise);

[ye1_cp, ye2_cp] = evea_design(ye_d,y_b_e1, y_b_e2, n_ofdm, n_cp, n_frame, n_L);
v_ye1_cp = mean(abs(ye1_cp));
v_ye2_cp = mean(abs(ye2_cp));