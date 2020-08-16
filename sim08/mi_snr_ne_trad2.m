function [v_y1_cp, v_y2_cp] = mi_snr_ne_trad2(data_ofdm, n_ofdm, n_cp ,n_frame, snr, cor)
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

h_12 = ray_model(d(3),taps(3));

h_12_t = ray_model_cor(h_12,cor,d(3));

pt = 10^(-2);
%snr = 100;
pw_noise = pt/10^(snr/10);

% all transmissions in ambient backscatter communication
y1_d = ofdm_trans(data_ofdm, h_12, pw_noise);
y1_d_t = ofdm_trans(data_ofdm, h_12_t, pw_noise);

[y1_cp] = receiver_design2(y1_d, n_ofdm, n_cp, n_frame, n_L);
[y2_cp] = receiver_design2(y1_d_t, n_ofdm, n_cp, n_frame, n_L);

v_y1_cp = mean(abs(y1_cp));
v_y2_cp = mean(abs(y2_cp));