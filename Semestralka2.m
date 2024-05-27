clc;
clear all;
close all;
M = 10000;
N = 100;
x = linspace(1, N);
x_1 = linspace(1, 95, 95);
x_cov = linspace(1, 6, 6);
%% prvni priklad
Q =3;
b = 0.5;
T = 1;
X = zeros(N, M);
W_odh = Q*(1-exp(-2*b*T));
for i=1:1:M
    X(1, i) = randn * sqrt(Q);
    for k=2:1:N
        X(k, i) = exp(-b*T)*X(k-1,i) + randn*sqrt(W_odh);
    end
end
figure;
hold on;
for n=1:1:5
    plot(x, X(:, n))
end
title('5 realizaci procesu')
xlabel('krok')
ylabel('hodnota X')
grid on;
hold off;
tau_max = 5;
cov_est = zeros(tau_max+1, 95);
for tau = 0:tau_max
    for k = 1:95
        sum_cov = 0;
        sum_cov = sum_cov + (X(k,:) - mean(X(k,:)))*(X(k+tau,:) - mean(X(k+tau,:)))';
        cov_est(tau+1, k) = sum_cov / M;
    end
end
cov_teor = zeros(6,95);
for j=1:1:6
    for i=1:1:95
        cov_teor(j,i) = exp(-b*T*(j-1))*Q;
    end
end
figure;
hold on;
grid on;
xlabel('krok')
ylabel('kovariacni funkce')
title('Namerena a teoriticka kovariacni funkce')
plot(x_1, cov_teor(1, :), 'black');
plot(x_1, cov_est(1, :));
plot(x_1, cov_teor(2, :), 'black', 'HandleVisibility', 'off');
plot(x_1, cov_est(2, :));
plot(x_1, cov_teor(3, :), 'black', 'HandleVisibility', 'off');
plot(x_1, cov_est(3, :));
plot(x_1, cov_teor(4, :), 'black', 'HandleVisibility', 'off');
plot(x_1, cov_est(4, :));
plot(x_1, cov_teor(5, :), 'black', 'HandleVisibility', 'off');
plot(x_1, cov_est(5, :));
plot(x_1, cov_teor(6, :), 'black', 'HandleVisibility', 'off');
plot(x_1, cov_est(6, :));
legend('teoriticka hodnota', 'tau 0','tau 1', 'tau 2', 'tau 3', 'tau 4', 'tau 5')
hold off;
%% druhy priklad
X = zeros(N, M);
W_odh = 1;
for i=1:1:M
    for k=2:1:N
        X(k, i)= X(k-1,i) + randn*sqrt(W_odh);
    end
end
figure;
hold on;
for n=1:1:5
    plot(x, X(:,n))
end
grid on;
xlabel('krok')
ylabel('Hodnota X')
title('5 realizaci procesu')
hold off;


Cov_est = zeros(tau_max+1, 95);
for tau = 0:tau_max
    for k = 1:95
        sum_cov = 0;
        sum_cov = sum_cov + (X(k,:) - mean(X(k,:)))*(X(k+tau,:) - mean(X(k+tau,:)))';
        Cov_est(tau+1, k) = sum_cov / M;
    end
end


cov_teor = zeros(1,95);
for j=1:1:95
    cov_teor(1,j) = j;
end

figure;
hold on;
grid on;
xlabel('krok')
ylabel('kovariacni funkce')
title('Namerena a teoriticka kovariacni funkce')
plot(x_1, cov_teor(1, :),'black');
plot(x_1, Cov_est(1, :));
plot(x_1, Cov_est(2, :));
plot(x_1, Cov_est(3, :));
plot(x_1, Cov_est(4, :));
plot(x_1, Cov_est(5, :));
plot(x_1, Cov_est(6, :));
legend('teoriticka hodnota', 'tau 0','tau 1', 'tau 2', 'tau 3', 'tau 4', 'tau 5')
hold off;
%% treti priklad
W_odch = 3;
V_odch = 2;
X0_mean = 1;
X0_odch = 5;
X = zeros(N, M);
Z = zeros(N, M);
for i=1:1:M
    X(1, i) = X0_mean + randn * sqrt(X0_odch);
    Z(1, i) = X(1,i) * 5 + randn*sqrt(V_odch);
    for k=2:1:N
        X(k, i) = 0.95*X(k-1, i) + 0.5*randn*sqrt(W_odch);
        Z(k, i) = 5*X(k, i) + randn*sqrt(V_odch);
    end
end
E_X = zeros(1, 100);
E_Z = zeros(1, 100);
for j=1:1:N
    E_X(1, j) = mean(X(j,:));
    E_Z(1, j) = mean(Z(j,:));
end
X_var = zeros(1,100);
Z_var = zeros(1,100);
for b=1:1:N
    sum_var_x = 0;
    sum_var_z = 0;
    for f=1:1:M
        sum_var_x = sum_var_x + (X(b, f) - E_X(1, b))^2;
        sum_var_z = sum_var_z + (Z(b, f) - E_Z(1, b))^2;
    end
    X_var(1,b) = sum_var_x/M;
    Z_var(1,b) = sum_var_z/M; 
end
E_teor = zeros(100, 2);
for k=1:1:100
    E_teor(k, 1) = 0.95^(k-1);
    E_teor(k, 2) = 5*0.95^(k-1);
end
VAR_teor = zeros(100, 2);
VAR_teor(1,1) = 5;
VAR_teor(1,2) = 5^3 + 2;
for k=2:1:100
    VAR_teor(k, 1) = (0.95^2)*VAR_teor(k-1, 1) + (0.5^2)*W_odch;
    VAR_teor(k, 2) = 25*VAR_teor(k, 1) + 2;
end


figure;
hold on;
plot(x, E_X(1, :));%Vykresleni stredni hodnoty x
plot(x, E_teor(:,1));
grid on;
title('Vykresleni stredni hodnoty vliciny x');
xlabel('krok');
ylabel('stredni hodnota x');
legend('namerena velicina', 'teoriticka velicina');
hold off;


figure;
hold on;
grid on;
title('Vykresleni stredni hodnoty vliciny z');
xlabel('krok');
ylabel('stredni hodnota z');
plot(x, E_Z(1, :));%Vykresleni streny hodnoty z
plot(x, E_teor(:, 2));
legend('namerena velicina', 'teoriticka velicina');
hold off;


figure;
hold on;
grid on;
title('Vykresleni odchylky vliciny x');
xlabel('krok');
ylabel('variace x');
plot(x, X_var(1, :));%Vykresleni odchylky x
plot(x, VAR_teor(:, 1));
legend('namerena hodnota', 'teorie');
hold off;


figure;
hold on;
grid on;
title('Vykresleni variace veliciny z');
xlabel('krok');
ylabel('variace z');
plot(x, Z_var(1, :));%Vykresleni odchylky z
plot(x, VAR_teor(:, 2));
legend('namerena hodnota', 'teorie');
hold off;

figure;
hold on;
title('Vykresleni realizace procesu pro velicinu Z');
for k=1:1:5
    plot(x, Z(:,k))
end
xlabel('krok')
ylabel('Hodnota Z')
grid on;
hold off;


figure;
hold on;
title('Vykresleni realizace procesu pro velicinu X');
for k=1:1:5
    plot(x, X(:,k))
end
xlabel('krok')
ylabel('Hodnota X')
grid on;
hold off;



