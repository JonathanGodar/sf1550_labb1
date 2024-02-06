%% 1a - plotta f(x)

f = @(x) x.^2 - 8 .* x - 12 .* sin(3.*x + 1) + 19;
f_prim = @(x) 2 .* x - 8 - 36 * cos(3 .* x + 1); 

x = [-10:0.1:10];

figure(1)
grid on
plot(x, f(x))

% Nollstallen:
% x = 2
% x = 2.7
% x = 4
% x = 4.8
% x = 6.2
% x = 6.7

% kod...

%% 1b - fixpunktiterationer

fixpunkt(2)
fixpunkt(2.7)
fixpunkt(4)
fixpunkt(4.8)
fixpunkt(6.2)
fixpunkt(6.7)

function rot = fixpunkt(gissning)
	tolerance = 1e-10;
    f_fixpunkt = @(x) 1/19 * (x^2 + 11 * x - 12 * sin(3 * x + 1)) + 1;
    
    x = gissning;
	x_prev = x + tolerance * 10;
    itterations = 0;
    while abs(x - x_prev) > tolerance
		x_prev = x;
        x = f_fixpunkt(x);
        itterations = itterations + 1;
    end
    
    disp("Gissning: ");
    gissning
    itterations
	rot = x
end

% kod...

%% 1c - Newton

% kod...

%% 1d - konvergensplottar

% figure(2);

% kod...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function xit = fixpunkt(x0,tau)

%  Indata:
%
%  x0  - startgissning (skalär)
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  xit - vektor av alla approximationer xit = [x0,x1,x2,x3,...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function xit = newton(x0,tau)

%  Indata:
%
%  x0  - startgissning (skalär)
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  xit - vektor av alla approximationer xit = [x0,x1,x2,x3,...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

% end
