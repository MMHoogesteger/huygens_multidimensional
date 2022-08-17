%% Mass inverse tester

syms MT Jp Jt
syms Ax Ay
syms Cx1 Cy1 Cp1
syms Cx2 Cy2 Cp2
syms Cx3 Cy3 Cp3

Mxy = diag([MT MT]);
Mt = diag([Jt,Jt,Jt]);

Axy = [Ax;Ay];
Cxy = [Cx1 Cx2 Cx3 ;Cy1 Cy2 Cy3];
Cp = [Cp1 Cp2 Cp3];

M = [Mxy  , Axy  , Cxy;
     Axy.', Jp   , Cp;
     Cxy.', Cp1.', Mt];
 
Mi = inv(M);
Mi = simplify(Mi,100);

