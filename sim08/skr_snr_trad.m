function  [v_y1_cp, v_y2_cp, v_ye1_cp, v_ye2_cp] = skr_snr_trad(data_ofdm, n_ofdm, n_cp ,n_frame, cor, d_e1, snr)

pt = 10^(-2);
%snr = 100;
alpha = 0.3 + 1i*0.4;
pw_noise = pt/10^(snr/10);

d = [8, 7, 3];
taps = [8, 7, 3];
n_L = max(taps(1),taps(2))+taps(3)-1;

h1 = ray_model(d(1),taps(1));
h2 = ray_model(d(2),taps(2));

h_12 = ray_model(d(3),taps(3));
h_21 = h_12;

% 假设h1信道和h3信道一摸一样，只考虑两个设备之间的信道，将其作为窃听信道
cor_h1_he = 1;  
cor_h12_he2 = cor;

% assume Eve locate closer to device 1
%cor = 0;
he = eave_ray_model(h1,cor_h1_he,d(1));
h_e2 = eave_ray_model(h_12,cor_h12_he2,d(3));

%d_e1 = 1;
tap_e1 = 1;
h_e1 = ray_model(d_e1,tap_e1);

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
