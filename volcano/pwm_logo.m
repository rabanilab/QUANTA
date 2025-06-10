function h = pwm_logo(PWM,N,AB,hdl)

if (nargin < 2)
    N = 1000;
end
if (nargin < 3)
    AB = 'ACGU';
end
if (nargin < 4)
    hdl = [];
end

% generate random sequences
rng('shuffle');
S = [];
for i = 1:size(PWM,2)
    c = cumsum(PWM(:,i));
    r = rand(N,1);
    p = zeros(size(r));
    for j = size(PWM,1):-1:1
        p(r<=c(j)) = j;
    end
    S = [S AB(p)'];
end

% display logo
if (isempty(hdl))
    [~,h] = SeqLogoFig(S,'CUTOFF',0.1);
else
    SeqLogoFig(S,'CUTOFF',0.1,'AXES_HANDLE',hdl);
end