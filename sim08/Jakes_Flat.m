function [h,tf]=Jakes_Flat(fd,Ts,Ns,t0,E0,phi_N)
% Inputs:
%       fd,Ts,Ns : Doppler frequency, sampling time, number of samples
%       t0, E0 : initial time, channel power
%       phi_N : inital phase of the maximum Doppler frequency sinusoid
% Outputs:
%       h, tf : complex fading vector, current time
if nargin<6
    phi_N=0; 
end
if nargin<5
    E0=1; 
end
if nargin<4
    t0=0; 
end
N0 = 8;         % As suggested by Jakes
N = 4*N0+2;     % an accurate approximation
wd = 2*pi*fd;   % Maximum Doppler frequency[rad]
t = t0+[0:Ns-1]*Ts; 
tf = t(end)+Ts; % Time vector and Final time
coswt=[sqrt(2)*cos(wd*t); 2*cos(wd*cos(2*pi/N*[1:N0]')*t)]; % Eq.(2.26)
h = E0/sqrt(2*N0+1) * exp(j*[phi_N pi/(N0+1)*[1:N0]])*coswt; % Eq.(2.23)

