function [y_cp] = receiver_design2(y, n_ofdm, n_cp, n_frame, n_L)

% calculate the desired information at receivering devices
y_cp = y(n_L:n_cp);
for n = 1:1:n_frame-1
    y_cp =[y_cp(:); (y(n*n_ofdm+n_L:n*n_ofdm+n_cp))];
end
