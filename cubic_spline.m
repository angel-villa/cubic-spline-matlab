% Angel Villa

clear; clc;

% rng(1)
x = [-1:0.02:2]';
y_actual = [sin(2*x.^2)];
y = y_actual + 0.2*randn(length(x),1);
N = length(x);

% Initialize lambda and EPE vectors with length N_l
test_lambda = 0:0.5:5;
N_l = length(test_lambda);
EPE = zeros(N_l,1);

% Calculate entries for EPE vector
for i=1:N_l
    EPE(i) = calc_E_val(x,y,test_lambda(i));
end

plot(test_lambda',EPE)
title('Plot of EPE at various values of lambda')
xlabel('lambda')
ylabel('EPE')

% Optimal lambda (lowest EPE) is 0
% large lambda gives a linear approximation
lambda = 0;

% Calculate knot points vector xi
xi = zeros(5,1);
for j=1:5
    index = floor(N/6*j);
    xi(j) = x(index);
end

% Basis matrix Nb
Nb = calc_Nb(N,x,xi);

% Calculate omega
omega = calc_omega(N,x,xi);

% Calculate theta vector
theta = (Nb'*Nb + lambda*eye(5)*omega)^-1*Nb'*y;

% Calculate f(X) for optimal lambda
f = zeros(N,1);
for i=1:N
    f(i) = Nb(i,:)*theta;
end

scatter(x,y,'b')
hold on
plot(x,f)
hold on
plot(x,y_actual,'k')
hold off
legend('y','f','actual function')
title('Plot of noisy data y and natural cubic spline spline f')
xlabel('x')

% Calculate EPE using leave one out cross-validation
function E = calc_E_val(x,y,lambda)
    E = 0;
    
    x_orig = x;
    y_orig = y;
    N_orig = length(x);
    
    N = length(x) - 1;
    
    % Iterate through 1:N, leaving out each index and testing
    % with that index each time
    for i=1:N_orig
        x = x_orig;
        y = y_orig;
        x(i) = [];
        y(i) = [];
        
        % Calculate knot points vector xi
        xi = zeros(5,1);
        for j=1:5
            index = floor(N/6*j);
            xi(j) = x(index);
        end
        
        % Basis matrix Nb
        Nb = calc_Nb(N,x,xi);
        
        % Calculate omega
        omega = calc_omega(N,x,xi);

        % Calculate theta vector
        theta = (Nb'*Nb + lambda*eye(5)*omega)^-1*Nb'*y;
        
        % Calculate bases N(x) for x(i) left out, for testing error
        Nt = calc_Nb(1,x_orig(i),xi);
        
        % Calculate testing error
        E = E + (y_orig(i) - Nt*theta)'*(y_orig(i) - Nt*theta);
    end
    
    E = E/N_orig;
end

% Calculate N basis matrix
function nb = calc_Nb(N,x,xi)
    nb = zeros(N,5);
    for i=1:N
        for j=1:5
            % Calculate dj(x) for each x(i), i=1:N
            d1 = (max((x(i)-xi(1))^3,0)-max((x(i)-xi(5))^3,0))/(xi(5)-xi(1));
            d2 = (max((x(i)-xi(2))^3,0)-max((x(i)-xi(5))^3,0))/(xi(5)-xi(2));
            d3 = (max((x(i)-xi(3))^3,0)-max((x(i)-xi(5))^3,0))/(xi(5)-xi(3));
            d4 = (max((x(i)-xi(4))^3,0)-max((x(i)-xi(5))^3,0))/(xi(5)-xi(4));
            if j==1
                nb(i,j) = 1;
            elseif j==2
                nb(i,j) = x(i);
            elseif j==3
                nb(i,j) = d1-d4;
            elseif j==4
                nb(i,j) = d2-d4;
            elseif j==5
                nb(i,j) = d3-d4;
            end
        end
    end
end
    
% Calculate N''i(x) to integrate for omega calculation
function Npp = calc_Npp(i,j,x,xi)
    Npp = 0;
    if j < 3
        Npp = 0;
    else
        Npp = 6*((max(x(i)-xi(j-2),0)-max(x(i)-xi(5),0)))/(xi(5)-xi(j-2));
        Npp = Npp - 6*((max(x(i)-xi(4),0)-max(x(i)-xi(5),0)))/(xi(5)-xi(4));
    end
end

% Calculate omega matrix
function omega = calc_omega(N,x,xi)
    omega = zeros(5,5);
    step = 0.005;
    for i=3:5
        for j=3:5
            for k=1:N
                Nppi = calc_Npp(k,i,x,xi);
                Nppj = calc_Npp(k,j,x,xi);
                omega(i,j) = omega(i,j) + Nppi*Nppj*step;
            end
        end
    end
end