function [v_y1_cp, v_y2_cp] = mi_snr_ne(data_ofdm, n_ofdm, n_cp ,n_frame, snr, cor)
% inputs:
%       data_ofdm: the modulated data with OFDM
%       n_ofdm:  the number of the ofdm symbols in a ofdm frame
%       n_cp:
%       n_frame: the number of the ofdm frames
%       snr: signal to noise rate (dB)
% outputs:
%       v_y1_cp: the ofdm symbols within CP part of devices 1
%       v_y2_cp: the ofdm symbols within CP part of devices 2

%% simulation parameters
d = [8, 7, 3];
taps = [8, 7, 3];
n_L = max(taps(1),taps(2))+taps(3)-1;

h1 = ray_model(d(1),taps(1));
h2 = ray_model(d(2),taps(2));
h_12 = ray_model(d(3),taps(3));

h1_t = ray_model_cor(h1,cor,d(1));
h2_t = ray_model_cor(h2,cor,d(2));
h_12_t = ray_model_cor(h_12,cor,d(3));


pt = 10^(-2);
%snr = 100;
alpha = 0.3 + 1i*0.4;
pw_noise = pt/10^(snr/10);

% all transmissions in ambient backscatter communication
y1_d = ofdm_trans(data_ofdm, h1, pw_noise);

y2_d = ofdm_trans(data_ofdm, h2, pw_noise);
back2 = ofdm_back(y2_d, alpha, n_ofdm, n_frame);
y1_b = ofdm_trans(back2, h_12, pw_noise);

y1_d_t = ofdm_trans(data_ofdm, h1_t, pw_noise);
back1_t = ofdm_back(y1_d_t, alpha, n_ofdm, n_frame);
y2_b_t = ofdm_trans(back1_t, h_12_t, pw_noise);

y2_d_t = ofdm_trans(data_ofdm, h2_t, pw_noise);


[y1,y1_cp] = receiver_design(y1_d, y1_b, n_ofdm, n_cp, n_frame, n_L);
[y2,y2_cp] = receiver_design(y2_d_t, y2_b_t, n_ofdm, n_cp, n_frame, n_L);

v_y1_cp = mean(abs(y1_cp));
v_y2_cp = mean(abs(y2_cp));