%% __Init__
format long
guessList = [2,2.7,4,4.8,6.2,6.7];
f = @(x) x.^2 - 8 .* x - 12 .* sin(3.*x + 1) + 19;
f_prim = @(x) 2 .* x - 8 - 36 * cos(3 .* x + 1); 

%% 1a - plotta f(x)
x = [-10:0.1:10];

figure(1)
plot(x, f(x))
grid on

% Nollstallen:
% x = 2
% x = 2.7
% x = 4
% x = 4.8
% x = 6.2
% x = 6.7

% kod...

%% 1b - fixpunktiterationer
for guess = guessList
    fixpunkt(guess)
end



% kod...

%% 1c - Newton
format long
for guess = guessList(1:5)
    output = newton(guess,f,f_prim,1e-10);
    display([guess,output{1},output{2}])

end
output = newton(guessList(6),f,f_prim,1e-10);
display([guessList(6),output{1},output{2}])
display(output{3})

%% 1d - konvergensplottar
% d1
xStarOut = newton(2,f,f_prim,1e-15);
xStar = xStarOut{1};

newtOut = newton(2,f,f_prim,1e-15);
fixOut = fixpunkt(2);

xNewt = newtOut{4};
xFix = fixOut{1};

errFix = abs(xFix-xStar);
errNewt = abs(xNewt-xStar);
xSize = 1:(size(errFix,2));
figure(2)
semilogy(xSize,[errNewt, zeros(size(errFix, 2) - size(errNewt, 2), 1)'],xSize,errFix)
grid on


%% d2
% Plotta e_n+1 som funktion av e_n
figure(3)
loglog(errNewt(1:end-1), errNewt(2:end), errFix(1:end-1), errFix(2:end))
grid on





% end
% b
function rot = fixpunkt(gissning)
	tolerance = 1e-10;
    f_fixpunkt = @(x) 1/19 * (x^2 + 11 * x - 12 * sin(3 * x + 1)) + 1;
    
    x = gissning;
	x_prev = x + tolerance * 10;
    itterations = 0;
    xList = [];
    while abs(x - x_prev) > tolerance
        xList = [xList,x];
		x_prev = x;
        x = f_fixpunkt(x);
        itterations = itterations + 1;
    end
    

    rot{1} = xList;
    rot{2} = itterations;

end
%%
% C
function output = newton(gissning,f,f_prim,tolerance)
    x = gissning;
    xPrev = x + tolerance * 2 ;
    list = [];
    numberOfIterations = 0;
    xList = [];
    while abs(x - xPrev) > tolerance
        xPrev = x;
        xList = [xList,x];
        x = x - f(x)/f_prim(x);
        numberOfIterations = numberOfIterations + 1;
        diff = abs(x-xPrev);
        list = [list,diff];
    end
    output{1}=x;
    output{2}=numberOfIterations;
    output{3}=list;
    output{4}=xList;
end
