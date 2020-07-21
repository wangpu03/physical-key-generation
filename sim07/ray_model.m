function h = ray_model(d,taps)
% generate the exponential rayleigh model
% inputs:
%       d: the distance between two objects
%       taps: the delay taps
% output:
%       h: the exponential Rayleigh model

lamda = 3;

h_cg = 5*d^(-lamda/2);

pow_h = exp(-(0:taps-1));
pow_h = h_cg*pow_h/norm(pow_h);

h = (randn(1,1) + 1i*randn(1,1)) * sqrt(pow_h);