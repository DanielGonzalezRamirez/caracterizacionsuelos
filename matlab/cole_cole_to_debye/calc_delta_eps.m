function [delta_eps] = calc_delta_eps(Nf, epsCC, sp, sf, eps_s, eps_inf, weights, rf)
%Calcula el t�rmino de offset para compensar la parte real de la funci�n Debye.
%   Nf: Tama�o del vector c. Corresponde al n�mero de muestras.
%   epsCC: Arreglo con los valores de la funci�n Cole-Cole
%   sp: Puntos de las muestras en frecuencia
%   sf: Frecuencias muestreadas
%   eps_s: Permitividad en frecuencia 0 de la funci�n Cole-Cole
%   eps_inf: Permitividad en frecuencia infinita de la funci�n Cole-Cole
%   weights: Pesos hallados para las funciones Debye
%   rf: Frecuencias de reflexi�n halladas para las funciones Debye

    samples = zeros(1,Nf);

    for i = 1:Nf
        poles = (eps_s - eps_inf) .* (weights./(1 + (sf(i)./rf).^2));
        samples(i) = real(epsCC(sp(i))) - eps_inf - sum(poles);
    end
    
    delta_eps = 1/Nf * sum(samples);
end

