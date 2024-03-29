%% __Init__
format long
guessList = [2,2.7,4,4.8,6.2,6.7];
f = @(x) x.^2 - 8 .* x - 12 .* sin(3.*x + 1) + 19;

%% 1a - plotta f(x)
x = [-10:0.1:14];

figure(1)
plot(x, f(x))
grid on

% Nollstallen från grafen:
% x = 2
% x = 2.7
% x = 4
% x = 4.8
% x = 6.2
% x = 6.7

%% 1b - fixpunktiterationer
tau = 1e-10;

disp("Fixpunkt med given funtion:")
guess = guessList(1);

xit = fixpunkt(guess, tau); 

fprintf("Hittade rot %f med startgissning: %f och Iterationer %d\n", xit(end), guess, size(xit, 2) - 1);
last_vals = abs(diff(xit));
disp(last_vals(end-9:end))


for guess = guessList(2:end)
    xit = fixpunkt(guess,tau);

    if abs(xit(end) - guess) < 1 
        fprintf("Hittade rot %f med startgissning: %f och Iterationer %d\n", xit(end), guess, size(xit, 2) - 1);
    end
end
    
%% 1c - Newton
tau = 1e-10;

disp("Newtons metod: ")
guess = guessList(1);

xit = newton(guess, tau);
fprintf("Hittade rot %f med startgissning: %f och Iterationer %d\n", xit(end), guess, size(xit, 2) - 1);
last_vals = abs(diff(xit));
disp(last_vals)


for guess = guessList(2:end)
    xit = newton(guess,tau);
    fprintf("Hittade rot %f med startgissning: %f och Iterationer %d\n", xit(end), guess, size(xit, 2) - 1);
    
end

%% 1d - konvergensplottar
% d1
x0 = 2;
tau = 1e-10;

x_newton15 = newton(x0,1e-15);
x_star = x_newton15(end);
x_fix = fixpunkt(x0, tau);

err_fix = abs(x_fix-x_star);
err_newt = abs(x_newton15-x_star);
x_vals_fix = 1:(size(err_fix,2));
x_vals_newt = 1:(size(err_newt, 2));

figure(2)
semilogy(x_vals_newt,err_newt,x_vals_fix,err_fix)
grid on

%% d2
% Plotta e_n+1 som funktion av e_n
figure(3)
loglog(err_newt(1:end-1), err_newt(2:end), err_fix(1:end-1), err_fix(2:end))
grid on
axis equal


%% b
function xit = fixpunkt(x0, tau)
    f_fixpunkt = @(x) 1/19 * (x^2 + 11 * x - 12 * sin(3 * x + 1)) + 1;
    x = x0;
	x_prev = x + tau * 10;
    xit = [x0];
    while abs(x - x_prev) >= tau 
		x_prev = x;
        x = f_fixpunkt(x);
        xit = [xit, x];
        
    end
end

%%
% C
function xit = newton(x0,tau)
    f = @(x) x.^2 - 8 .* x - 12 .* sin(3.*x + 1) + 19;
    f_prim = @(x) 2 .* x - 8 - 36 * cos(3 .* x + 1); 

    x = x0;
    % Skulle kunna utskiftas om x0 stor och tolerance är liten. Alternativ är att använda inf.
    xPrev = x + tau * 10 ;
    xit = [xPrev];
    while abs(x - xPrev) > tau
        xPrev = x;
        x = x - f(x)/f_prim(x);
        xit = [xit,x];
        
        if size(xit, 2) > 4000
            break
        end
    end
end
