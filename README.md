%a) 

clc, clearvars

C_grid = linspace(0,200,20000);

%b) 

sigma_grid = linspace(1.1,50,1000);

%c)

function [Utilidad] = Utilidad(c, beta, sigma)
    c1 = c(1);
    c2= c(2);
    
    if c1 >0 && c2 > 0
        Utilidad = ((c1.^(1 - sigma) - 1) ./ (1 - sigma)) + beta * ((c2.^(1 - sigma) - 1) ./ (1 - sigma));
    else
        Utilidad = -1e10;
    end
end

%%
%d) 
tic
beta = 0.9;

R = 1.03;

y_total = 100;

consumo = zeros(1000, 20000); 

C = zeros(20000,2);

for i = 1:1000
    sigma = sigma_grid(i);
    for j = 1:20000

        c1 = C_grid(j);
        c2 = (y_total - c1) * (R);

        c = [c1, c2];
        consumo(i, j) = Utilidad(c, beta, sigma);
    end 
end 
toc

%%
%e)
tic
maximo = zeros(1,1000);
posicion = zeros(1,1000);
consumo_optimo = zeros(1000,2); % Cambio de tamaño (1000 filas, 2 columnas)

for i = 1:1000
    [maximo(i), posicion(i)] = max(consumo(i,:)); 

    c1_optimo = C_grid(posicion(i));
    
    consumo_optimo(i,1) = c1_optimo;
    consumo_optimo(i,2) = (y_total - c1_optimo) * R; 
end
toc

%%
%e) 
tic
perfil_consumo = zeros(1000,1);

for i = 1:1000
    perfil_consumo(i) = consumo_optimo(i,2)/consumo_optimo(i,1);
end

figure(1)

plot(sigma_grid, perfil_consumo)
grid on
xlabel("sigma") 
ylabel("Consumo Relativo (c2/c1)") 
title("Perfil de Consumo Óptimo para cada valor de Sigma")
toc
