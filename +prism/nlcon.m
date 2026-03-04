function [c, ceq] = nlcon(mixing, mixingGuess)

r0 = mixingGuess(:,1)./mixingGuess(:,2);
r = mixing(:,1)./mixing(:,2);

cd = -r0./r +0.5;
cu = r0./r -2;

c = cat(1, cd, cu);
ceq = [];
end