% =========================================================================
% Ejercicio 1 - Clase de Laboratorio: Integración Numérica
% =========================================================================
clc; clear; close all;

disp('--- RESOLUTOR DE INTEGRALES NUMÉRICAS ---');
disp('Ingresa la función de esta forma:');
disp('1) 1/sqrt(2*pi)*exp(-x^2/2)');
disp('2) cos(x^2)');
disp('3) x^5 - 2*x^3 + 4');
disp('-----------------------------------------');

% 1. Ingreso de datos
str_f = input('Ingrese la función f(x): ', 's');
a = input('Ingrese el límite inferior a: ');
b = input('Ingrese el límite superior b: ');

syms x;
f_sym = str2sym(str_f);
f_func = matlabFunction(f_sym);

% 2. Identificación de polinomio y ajuste de 'n' para Gauss
n_gauss = 5; % Valor por defecto según la diapositiva
es_polinomio = false;

try
    coeficientes = sym2poly(f_sym); 
    grado = length(coeficientes) - 1;
    % Fórmula para solución exacta: 2n - 1 >= grado -> n >= (grado + 1)/2
    n_gauss = ceil((grado + 1) / 2);
    es_polinomio = true;
    fprintf('\n[!] Polinomio de grado %d detectado. Gauss usará n = %d para ser exacto.\n\n', grado, n_gauss);
catch
    fprintf('\n[!] La función no es polinómica. Gauss usará n = %d nodos.\n\n', n_gauss);
end

% Integral exacta para sacar el error absoluto
I_exacta = double(int(f_sym, x, a, b));

% 3. Cálculos con los métodos (N = 12)
N = 12;

[I_trap, ev_trap] = trapecio_compuesto(f_func, a, b, N);
[I_s13, ev_s13]   = simpson_13_compuesto(f_func, a, b, N);
[I_s38, ev_s38]   = simpson_38_compuesto(f_func, a, b, N);
[I_gl, ev_gl]     = gauss_legendre(f_func, a, b, n_gauss);

% Errores absolutos
E_trap = abs(I_exacta - I_trap);
E_s13  = abs(I_exacta - I_s13);
E_s38  = abs(I_exacta - I_s38);
E_gl   = abs(I_exacta - I_gl);

% 4. Mostrar resultados en tabla
disp('Los resultados deben mostrarse así:');
fprintf('Método                | Puntos F(x) | Aproximación   | Error Absoluto\n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('Trapecio Compuesto    | %11d | %14.8f | %14.4e\n', ev_trap, I_trap, E_trap);
fprintf('Simpson 1/3 Compuesto | %11d | %14.8f | %14.4e\n', ev_s13, I_s13, E_s13);
fprintf('Simpson 3/8 Compuesto | %11d | %14.8f | %14.4e\n', ev_s38, I_s38, E_s38);
fprintf('Gauss-Legendre        | %11d | %14.8f | %14.4e\n', ev_gl, I_gl, E_gl);


% =========================================================================
% FUNCIONES DE INTEGRACIÓN
% =========================================================================

function [I, evals] = trapecio_compuesto(f, a, b, N)
    h = (b-a)/N;
    X = a:h:b;
    Y = f(X);
    I = (h/2) * (Y(1) + 2*sum(Y(2:end-1)) + Y(end));
    evals = N + 1; % 13 puntos para N=12
end

function [I, evals] = simpson_13_compuesto(f, a, b, N)
    h = (b-a)/N;
    X = a:h:b;
    Y = f(X);
    I = (h/3) * (Y(1) + 4*sum(Y(2:2:end-1)) + 2*sum(Y(3:2:end-2)) + Y(end));
    evals = N + 1;
end

function [I, evals] = simpson_38_compuesto(f, a, b, N)
    h = (b-a)/N;
    X = a:h:b;
    Y = f(X);
    I = (3*h/8) * (Y(1) + 3*sum(Y(2:3:end-2)) + 3*sum(Y(3:3:end-1)) + 2*sum(Y(4:3:end-3)) + Y(end));
    evals = N + 1;
end

function [I, evals] = gauss_legendre(f, a, b, n)
    % Cálculo algorítmico de raíces y pesos (Matriz tridiagonal)
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T); % Eigenvalores = raíces, Eigenvectores = pesos
    
    raices = diag(D);
    [raices, idx] = sort(raices);
    pesos = 2 * V(1, idx).^2;
    pesos = pesos(:);
    
    % Cambio de variable a los límites [a, b]
    x_map = ((b-a)*raices + (a+b))/2;
    Y = f(x_map);
    
    I = ((b-a)/2) * sum(pesos .* Y);
    evals = n; % Solo evalúa 'n' veces
end