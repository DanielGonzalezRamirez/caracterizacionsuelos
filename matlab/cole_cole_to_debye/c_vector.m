function [c] = c_vector(Nf, epsCC, sp, eps_s, eps_inf)
%Calcula el vector c para la optimización lineal de los pesos de las funciones Debye.
%   Nf: Tamaño del vector c. Corresponde al número de muestras.
%   epsCC: Arreglo con los valores de la función Cole-Cole
%   sp: Puntos de las muestras en frecuencia
%   eps_s: Permitividad en frecuencia 0 de la función Cole-Cole
%   eps_inf: Permitividad en frecuencia infinita de la función Cole-Cole
    c = zeros(Nf,1);
    for i = 1:Nf
        c(i) = abs(imag(epsCC(sp(i))))/(eps_s - eps_inf);
    end
end

