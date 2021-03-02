function [D] = D_matrix(Nf, Np, sf, rf)
%Calcula la matriz D para la optimización lineal de los pesos de las funciones Debye.
%   Nf: # de filas de la matriz D. Corresponde al número de muestras.
%   Np: # de columnas de la matriz D. Corresponde al número de polos.
%   sf: Vector con frecuencias de muestreo.
%   rf: Vector con frecuencias de relajación.
    D = zeros(Nf, Np);
    for i = 1:Nf
        for j = 1:Np
            D(i,j) = (sf(i)/rf(j)) / (1+(sf(i)/rf(j))^2);
        end
    end
end

