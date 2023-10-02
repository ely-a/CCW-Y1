function [M] = getM(m, rho, c, xc, xf, Ialpha) 

%function which returns the mass matrix (M) of the aeroelastic

%use help function for getK() to see inputs
%outputs: M      -   mass matrix of dimensions 8x8
%inputs:  m      -   mass of the airfoil per unit depth
%         Ialpha -   sectional polar second moment of intertia per unit
%                    depth (m)
%         c      -   airfoil chord length (m)
%         xc     -   sectional cg position from l.e. (% chord)
%         xf     -   sectional flexural centre position from l.e. (% chord)
%         rho    -   density of air (kg/m^3)

%creating a general matrix of 0
M = eye(8);

%defining M11, M12, M21 and M22
M(1,1) = -m-rho*pi*(c/2)^2;
M(1,2) = rho*pi*(xf*c - c/2)*(c/2)^2 - (xc - xf)*c*m;
M(2,1) = M(1,2);
M(2,2) = - Ialpha - rho*pi*((c/2)^2)*(xf*c - c/2)^2 + (1/8)*(c/2)^2;

end 
