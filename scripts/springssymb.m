%% Derive springs etc

syms W L x y phi syms Ls Fz
assume(L >0)
assume(W >0)

rur = [x + L*cos(phi) - W*sin(phi); y + L*sin(phi) + W*cos(phi)];
rul = subs(rur,L,-L);
rlr = subs(rur,W,-W);
rll = subs(rlr,L,-L);

rurf = [L;W];
rulf = [-L;W];
rlrf = [L;-W];
rllf = [-L;-W];

cm = [x;y];

aur = rurf-rur;
aul = rulf-rul;
alr = rlrf-rlr;
all = rllf-rll;

magaur = norm(aur);
magaul = norm(aul);
magalr = norm(alr);
magall = norm(all);

thetaur = asin(magaur/Ls);
thetaul = asin(magaul/Ls);
thetalr = asin(magalr/Ls);
thetall = asin(magall/Ls);

Fur = tan(thetaur)*Fz;
Ful = tan(thetaul)*Fz;
Flr = tan(thetalr)*Fz;
Fll = tan(thetall)*Fz;

Furvec = Fur*(aur/magaur);
Fulvec = Ful*(aul/magaul);
Flrvec = Flr*(alr/magalr);
Fllvec = Fll*(all/magall);

Frt = Furvec+Fulvec+Flrvec+Fllvec;


lur = rur-cm;
lul = rul-cm;
llr = rlr-cm;
lll = rll-cm;

Mur = -Furvec(2)*lur(1)+Furvec(1)*lur(2);
Mul = -Fulvec(2)*lul(1)+Fulvec(1)*lul(2);
Mlr = -Flrvec(2)*llr(1)+Flrvec(1)*llr(2);
Mll = -Fllvec(2)*lll(1)+Fllvec(1)*lll(2);

Mt = Mur+Mul+Mlr+Mll;
%%
Frt = simplify(Frt,100);
Mt = simplify(Mt,100);

%%
xs = 0;
ys = 0;
phis = 0;

Frtxdx = diff(Frt(1),x);
Frtxdy = diff(Frt(1),y);
Frtxdphi = diff(Frt(1),phi);

Frtxdxs = subs(Frtxdx,{x,y,phi},{xs,ys,phis})
Frtxdys = subs(Frtxdy,{x,y,phi},{xs,ys,phis})
Frtxdphis = subs(Frtxdphi,{x,y,phi},{xs,ys,phis})

%
Frtydx = diff(Frt(2),x);
Frtydy = diff(Frt(2),y);
Frtydphi = diff(Frt(2),phi);

Frtydxs = subs(Frtydx,{x,y,phi},{xs,ys,phis})
Frtydys = subs(Frtydy,{x,y,phi},{xs,ys,phis})
Frtydphis = subs(Frtydphi,{x,y,phi},{xs,ys,phis})

%
Mtdx = diff(Mt,x);
Mtdy = diff(Mt,y);
Mtdphi = diff(Mt,phi);

Mtdxs = subs(Mtdx,{x,y,phi},{xs,ys,phis})
Mtdys = subs(Mtdy,{x,y,phi},{xs,ys,phis})
Mtdphis = subs(Mtdphi,{x,y,phi},{xs,ys,phis})