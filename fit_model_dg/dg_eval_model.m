function x = dg_eval_model(t,P)
% P = <logX0,dg,t0,t1>

x = zeros(size(t));
i = t<=P(3);
x(i) = P(1);
i = (t>P(3)).*(t<P(4))==1;
x(i) = P(1) - P(2)*(t(i)-P(3));
i = t>=P(4);
x(i) = P(1) - P(2)*(P(4)-P(3));
