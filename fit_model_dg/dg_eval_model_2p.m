function x = dg_eval_model_2p(t,P)
% P = <logX0,dg1,dg2,t0,t1>

x = zeros(size(t));
i = t<=P(4);
x(i) = P(1) - P(2)*t(i);
i = (t>P(4)).*(t<P(5))==1;
x(i) = (P(1) - P(2)*P(4)) - P(3)*(t(i)-P(4));
i = t>=P(5);
x(i) = (P(1) - P(2)*P(4)) - P(3)*(P(5)-P(4));
