function [y,y_cp] = receiver_design(y_d, y_b, n_ofdm, n_cp, n_frame, n_L)
% combinaiton of two signals

sym_rem = mod(n_ofdm-mod(length(y_d),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
y_d_padded = [y_d;padding];

sym_rem = mod(n_ofdm-mod(length(y_b),n_ofdm),n_ofdm);
padding = repmat(0+0i,sym_rem,1);
y_b_padded = [y_b;padding];

y = y_d_padded+y_b_padded;

% calculate the desired information at receivering devices
n_fft = n_ofdm - n_cp;
zb = y(n_L:n_cp) - y(n_fft+n_L:n_fft+n_cp);
for n = 1:1:n_frame-1
    zb =[zb(:); (y(n*n_ofdm+n_L:n*n_ofdm+n_cp)-y(n*n_ofdm+n_fft+n_L:n*n_ofdm+n_fft+n_cp))];
end

zd = y(n_L:n_cp) + y(n_fft+n_L:n_fft+n_cp);
for n = 1:1:n_frame-1
    zd =[zd(:); (y(n*n_ofdm+n_L:n*n_ofdm+n_cp)+y(n*n_ofdm+n_fft+n_L:n*n_ofdm+n_fft+n_cp))];
end
y_cp = zb .* zd;
