function drv = CRTBP(t, rv)
% ------------------------------------------------------------------------
% Inputs 
%   t = [Nx1] time vector (orbit is Keplerian, doesn't matter) 
%   x = [6x1] state vector 
% 
% Outputs 
%   dx = [6x1] derivative of state vector 
% ------------------------------------------------------------------------

global mu 

drv = zeros(6, 1);   % force column vector 

x   = rv(1); 
y   = rv(2); 
z   = rv(3); 
dx  = rv(4); 
dy  = rv(5); 
dz  = rv(6); 

% r_norm  = sqrt(x^2 + y^2 + z^2); 
r1  = sqrt( (x+mu)^2 + y^2 + z^2 ); 
r2  = sqrt( (x-1+mu)^2 + y^2 + z^2 ); 

C1  = (1-mu)/r1^3; 
C2  = mu/r2^3; 

ddx = 2*dy + x - C1*(x+mu) - C2*(x-1+mu); 
ddy = -2*dx + y - C1*y - C2*y; 
ddz = -C1*z - C2*z; 

drv(1:6) = [dx; dy; dz; ddx; ddy; ddz]; 

end