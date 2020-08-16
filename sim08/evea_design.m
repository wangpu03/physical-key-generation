function [ye1_cp, ye2_cp] = evea_design(y_d,y_b_e1, y_b_e2, n_ofdm, n_cp, n_frame, n_L)
%

sym_rem = mod(n_ofdm-mod(length(y_d),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
y_d_padded = [y_d;padding];

sym_rem = mod(n_ofdm-mod(length(y_b_e1),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
y_b_e1_padded = [y_b_e1;padding];

sym_rem = mod(n_ofdm-mod(length(y_b_e2),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
y_b_e2_padded = [y_b_e2;padding];

y_e2 = y_d_padded+y_b_e2_padded;

y_e1 = y_d_padded+y_b_e1_padded;
% calculate the desired information at receivering devices
n_fft = n_ofdm - n_cp;
zb_e2 = y_e2(n_L:n_cp) - y_e2(n_fft+n_L:n_fft+n_cp);
for n = 1:1:n_frame-1
    zb_e2 =[zb_e2(:); (y_e2(n*n_ofdm+n_L:n*n_ofdm+n_cp)-y_e2(n*n_ofdm+n_fft+n_L:n*n_ofdm+n_fft+n_cp))];
end

zd_e2 = y_e2(n_L:n_cp) + y_e2(n_fft+n_L:n_fft+n_cp);
for n = 1:1:n_frame-1
    zd_e2 =[zd_e2(:); (y_e2(n*n_ofdm+n_L:n*n_ofdm+n_cp)+y_e2(n*n_ofdm+n_fft+n_L:n*n_ofdm+n_fft+n_cp))];
end
ye2_cp = zb_e2 .* zd_e2;  %%发射信号和直射信号相乘，得到第一类窃听方法

n_fft = n_ofdm - n_cp;
zb_e1 = y_e1(n_L:n_cp) - y_e1(n_fft+n_L:n_fft+n_cp);
for n = 1:1:n_frame-1
    zb_e1 =[zb_e1(:); (y_e1(n*n_ofdm+n_L:n*n_ofdm+n_cp)-y_e1(n*n_ofdm+n_fft+n_L:n*n_ofdm+n_fft+n_cp))];
end

ye1_cp = zb_e1 .* zb_e2;  %%两个反射信号相乘，得到第二类窃听方法
