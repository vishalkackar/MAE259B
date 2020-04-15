function F = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, l_k, EI)
%
% This function returns the derivative of bending energy E_k^b with respect
% to x_{k-1}, y_{k-1}, x_k, y_k, x_{k+1}, and y_{k+1}.
%
F = zeros(6,1);

F(1) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((ykp1 - yk) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (-xkp1 + xk)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);
F(2) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((-xkp1 + xk) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (-ykp1 + yk)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);
F(3) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((-ykp1 + ykm1) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (xkp1 - 0.2e1 * xk + xkm1)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);
F(4) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((-xkm1 + xkp1) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (ykp1 - 0.2e1 * yk + ykm1)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);
F(5) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((yk - ykm1) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (xk - xkm1)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);
F(6) = 0.4e1 * tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) * ((-xk + xkm1) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) - (-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 * (yk - ykm1)) / ((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) ^ 2 / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk)) ^ 2 + 0.1e1) * (0.1e1 + tan(atan2((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) , ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1) ^ 2);

F = 0.5 * EI * F / l_k;

end