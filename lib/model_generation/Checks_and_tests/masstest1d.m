N = 10;
syms MT MR MR2

theta = sym('theta',[N,1]);
%dtheta = sym('dtheta',[N,1]);

Mh = sym(zeros(1+N,1+N));

Mh(1,1) = MT;

for met_k = 1:N
    Mh(1,1+met_k) = MR*cos(theta(met_k));
    Mh(1+met_k,1) = Mh(1,1+met_k);
    
    Mh(1+met_k,1+met_k) = MR2;
end

MM = Mh;

spy(inv(MM))