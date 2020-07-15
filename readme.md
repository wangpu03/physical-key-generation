this document is used to discuss the mutual information between two variables

# 1. mi_variables
% calculate the mutual information between two variables X, Y, that multiply from two Guassian variables
% f channel information (Guassin variables)
% g channel information (Guassin variables)
% z1 and z1 noise signals ((Guassin variables))
% X = f.*g + z1;
% Y = f.*g + z2;

# 2. mi_of_norm
% calculate the mutual information between two Guassin variables x, y
% r_h channel information (Guassin variables)
% z1 and z1 noise signals ((Guassin variables))
% x = r_h+z1;
% y = r_h+z2;

# 3. entropy_of_norm
% calculate the entropy of the norm variable with different manners.

# mi_BC_system_h12_norm.m
% based on the direct channel model to calculate the mutual information between two backscatter devices
% calculate the theortical and real mutual information between two
% variables (v12,v21) from ambient backscatter system model
% h12,f1,f2: three channel variables of the corresponding channel information(Guassin variables)
% z1,z2: the noise information(Guassin variables)
% v12 = 4*alpha*(rho^2)*Pt*h12*f1.*f2+z1;
% v21 = 4*alpha*(rho^2)*Pt*h12*f2.*f1+z2;

# mi_BC_system_h12_cons.m
% calculate the theortical and real mutual information between two
% variables (v12,v21) from ambient backscatter system model
% f1,f2: three channel variables of the corresponding channel information(Guassin variable)
% h12 the channe information between two device (scale variable)
% z1,z2: the noise information(Guassin variable)
% v12 = 4*alpha*(rho^2)*Pt*h12*f1.*f2+z1;
% v21 = 4*alpha*(rho^2)*Pt*h12*f2.*f1+z2;

# key-gen-sim-01 and key-gen-sim-02
% simulate the key generation scheme with the based channel model that are modeled as a Gaussian variables. These two documents are mainly used to verify the idea and mutual information between two devices.


