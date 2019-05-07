function [c_real,c_theo,c_sk, energy] = key_gen_func(rho, Pt)

%  set the parameters
%  channels, backscatter coefficient, energy harvesting coefficient and
%  the transmit power
%  the form of the channel model f=g*d^(L/2)--(d-distance,g-CSCG, f-channel)

alpha = 0.8;        %the constant of backscatter coefficient
%Pt = pt;           %transmission power of the RF source
%rho = 0.2;          %power split ratio for information receiver, range is [0,1]

sigma_f1 = sqrt(0.868);    %channel f1
f1 = normrnd(0,sigma_f1,[1 1000000])';

sigma_f2 = sqrt(0.680);    %channel f2
f2 = normrnd(0,sigma_f2,[1 1000000])';

sigma_h12 = sqrt(1.93);          %channel h12 or h21
h12 = normrnd(0,sigma_h12,[1 1000000])';

sigma_z1 = 1;       % noise at the device BD1
z1 = normrnd(0,sigma_z1,[1 1000000])';

sigma_z2 = 1;       % noise at the device BD2
z2 = normrnd(0,sigma_z2,[1 1000000])';

v1 = 4*alpha*(rho^2)*Pt*h12.*f1.*f2;
v2 = 4*alpha*(rho^2)*Pt*h12.*f2.*f1;
v12 = v1+z1;
v21 = v2+z2;
%v12 = 4*alpha*(rho^2)*Pt*h12*f1.*f2+z1;
%v21 = 4*alpha*(rho^2)*Pt*h12*f2.*f1+z2;

c_theo = 1/2*log2(1+(var(v1)/(var(z1)+var(z2)+(var(z1)*var(z2)/var(v1)))));
c_real = mi(v12,v21);

sigma_h_1e = sqrt(5.4);
h_1e = normrnd(0,sigma_h_1e,[1 1000000])';
sigma_h_2e = sqrt(5.4);
h_2e = normrnd(0,sigma_h_2e,[1 1000000])';

v_e = 4*alpha^2*Pt*h_1e.*h_2e.*f1.*f2;
c_sk = cmi(v12,v21,v_e);

zeta = 0.8;
energy  = 1/2*zeta*(1-rho)^2*(var(f1)+var(f2))*Pt;