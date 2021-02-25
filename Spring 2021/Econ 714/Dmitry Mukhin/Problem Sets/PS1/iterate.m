function [Kp,Cp] = iterate(K, C, D, psigma, palpha, pbeta, pdelta)
Kp = K.^palpha + (1-pdelta)*K - C - D;
Cp = C.*(pbeta * (palpha*K.^(palpha - 1) + 1 - pdelta)).^(1/psigma);
end