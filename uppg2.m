%punkter([0.1,0.9,2.1,0.9]',1,1,[0,0],[2,0])
X = punkter([0.1,4,2.1,0.9]',1.5,0.8,[-1.5,3],[1,1])
hold on
fplot(@(t) 1.5*cos(t)-1.5, @(t) 1.5*sin(t)+3)
hold on
axis equal
fplot(@(t) 0.8*cos(t)+1, @(t) 0.8*sin(t)+1)
hold on
line([X(1) X(3)],[X(2) X(4)])
% Kod för uppgift 2... 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,X]=langd(ra,rb,rc,a,b,c)

%  Indata:
%
%  ra - radie för cirkel a (skalär)
%  rb - radie för cirkel a (skalär)
%  rc - radie för cirkel a (skalär)
%  a  - kolumnvektor med koordinater för mittpunkt a (xa,ya)^T
%  b  - kolumnvektor med koordinater för mittpunkt b (xb,yb)^T
%  c  - kolumnvektor med koordinater för mittpunkt c (xc,yc)^T
%
%  Utdata:
%
%  L - längden på snöret (skalär)
%  X - De tre lösningsvektorerna samlade i en 
%      matris (4 rader, 3 kolumner), X=[Xab, Xbc, Xca]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

end
%%
function [Xrot, iter] = punkter(X0,ra,rb,a,b)
    a1 = a(1);
    a2 = a(2);
    b1 = b(1);
    b2 = b(2);

		% x1 <- x1(1)
		% y1 <- x1(2)

		% x2 <- x2(1)
		% y2 <- x2(2)

    F = @(x) [(x(1)-a1)^2+(x(2)-a2)^2-ra^2;
              (x(3)-b1)^2+(x(4)-b2)^2-rb^2;
              (x(1)-x(3))*(x(1)-a1)+(x(2)-x(4))*(x(2)-a2);
              (x(1)-x(3))*(x(3)-b1)+(x(2)-x(4))*(x(4)-b2)];

    jacobian = @(x) [2*(x(1)-a1) 2*(x(2)-a2) 0 0;
                     0 0 2*(x(3)-b1) 2*(x(4) -b2);
                     2*x(1)-a1-x(3) 2*x(2)-a2-x(4) a1-x(1) a2-x(2);
                     x(3)-b1 x(4)-b2 x(1)+b1-2*x(3) x(2)+b2-2*x(4)];

		x = X0;

        xList = [x]
		tolerance = 1e-10;
        iterations = 0;
		xPrev = x + 2 * tolerance;
		while max(abs(x - xPrev)) > tolerance
			xPrev = x;
			% invJacobian = inv(jacobian(x));
			x = x - jacobian(x) \ F(x);
            iterations = iterations + 1;
            xList = [xList, x]
        end
        disp(x);
        disp(iterations);
        Xrot = x
        iter = iterations
        plot(xList(1:4:iterations * 4), ...
            xList(2:4:iterations * 4), ...
            xList(3:4:iterations * 4), ...
            xList(4:4:iterations * 4), '+')
end

    % jacobian = @(x1, x2) [2*(x1(1)-a(1)) 2*(x1(2)-a2) 0 0;
    %                            0 0 2*(x2(1)-b(1)) 2*(x2(2) -b2);
    %                            2*x1(1)-a(1)-x2(1) 2*x1(2)-a2-x2(2) a(1)-x1(1) a2-x1(2);
    %                            x2(1)-b(1) x2(2)-b2 x1(1)+b(1)-2*x2(1) x1(2)+b2-2*x2(2)];

    % jacobian_fungerande = @(x1,y1,x2,y2) [2*(x1-a1) 2*(y1-a2) 0 0;
    %                            0 0 2*(x2-b1) 2*(y2 -b2);
    %                            2*x1-a1-x2 2*y1-a2-y2 a1-x1 a2-y1;
    %                            x2-b1 y2-b2 x1+b1-2*x2 y1+b2-2*y2];
    %  Indata:
%
%  X0 - kolumnvektor med startgissningarna (x1,y1,x2,y2)^T
%  ra - radie för cirkel a (skalär)
%  rb - radie för cirkel a (skalär)
%  a  - kolumnvektor med koordinater för mittpunkt a (xa,ya)^T
%  b  - kolumnvektor med koordinater för mittpunkt b (xb,yb)^T
%
%  Utdata:
%
%  Xrot - kolumnvektor med lösningen (x1,y1,x2,y2)^T
%  iter - antal iterationer som använts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...