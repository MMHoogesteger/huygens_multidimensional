function [U] = escapement_ham(mdl,dtheta,Htheta,HStar,U_gains)
U = sym(zeros(mdl.Nq,1));

U(mdl.NqP+1:mdl.NqP+mdl.N) = U_gains.*dtheta.*(Htheta-HStar);

end