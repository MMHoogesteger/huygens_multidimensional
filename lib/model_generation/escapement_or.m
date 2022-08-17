function [U] = escapement_or(mdl,theta, dtheta, epsilon, theta_s, theta_e, tau)
U = sym(zeros(mdl.Nq,1));

si = tanh(dtheta./epsilon);
U(mdl.NqP+1:mdl.NqP+mdl.N) = tau.*si.*(tanh((si.*theta-theta_s)./epsilon)-tanh((si.*theta-theta_e)./epsilon))/2;

end