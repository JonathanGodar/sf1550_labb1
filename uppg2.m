%punkter([0.1,0.9,2.1,0.9]',1,1,[0,0],[2,0])
aCoords = [-1; 1.5];
bCoords = [3; 0.5];
cCoords = [0; -2];

aRadius = 1;
bRadius = 1.2;
cRadius = 1.7;

tic;
a_b_guess = aCoords + [0;1];
b_a_guess = bCoords + [0;1];

a_b_coords= punkter([a_b_guess; b_a_guess], aRadius, bRadius, aCoords, bCoords);

a_c_guess = aCoords + [-1;0];
c_a_guess = cCoords + [-1;0];
a_c_coords = punkter([a_c_guess; c_a_guess], aRadius, cRadius, aCoords, cCoords);

b_c_guess = bCoords + [0;-1];
c_b_guess = cCoords + [1;0];
b_c_coords = punkter([b_c_guess; c_b_guess], bRadius, cRadius, bCoords, cCoords);
toc

a_to_b = a_b_coords(1:2);
b_to_a = a_b_coords(3:4);

a_to_c = a_c_coords(1:2);
c_to_a = a_c_coords(3:4);

b_to_c = b_c_coords(1:2);
c_to_b = b_c_coords(3:4);


% line([a_to_b(1) b_to_a(1)], [a_to_b(2) b_to_a(2)])



drawSphere(aCoords, aRadius, 'b');
hold on 
axis equal
drawSphere(bCoords, bRadius, 'b');
drawSphere(cCoords, cRadius, 'b');

plot_line_from_matrix([a_to_b, b_to_a]);
plot_line_from_matrix([a_to_c, c_to_a]);
plot_line_from_matrix([c_to_b, b_to_c]);
hold off


points_list = [a_to_b, a_to_c, b_to_a, b_to_c, c_to_a, c_to_b];
circle_center_list = [aCoords, bCoords, cCoords];
radius_list = [aRadius, bRadius, cRadius];

disp("Svar: ")
calculate_rope_length(circle_center_list, radius_list, points_list)

% a_perim = calculate_perimiter(aCoords, aRadius, b_to_c, b_to_a);
% b_perim = calculate_perimiter(bCoords, bRadius, b_to_c, b_to_a);

% line([a_to_b(1) a_to_b(3)],[a_to_b(2) a_to_b(4)]);
% line([a_to_c(1) a_to_c(3)],[a_to_c(2) a_to_c(4)]);
% line([b_to_c(1) b_to_c(3)],[b_to_c(2) b_to_c(4)]);
% hold off

% calculate_perimiter(bCoords, bRadius, b_to_c, b_to_a)


% b_to_a

% a_to_c
% c_to_a

% b_to_c
% c_to_b


% drawSphere(bCoords, bRadius, 'b')
% drawSphere(cCoords, cRadius,'g')



% X = punkter([0.1,4,2.1,0.9]',1.5,0.8,[-1.5,3],[1,1])
% hold on
% fplot(@(t) 1.5*cos(t)-1.5, @(t) 1.5*sin(t)+3)
% hold on
% axis equal
% fplot(@(t) 0.8*cos(t)+1, @(t) 0.8*sin(t)+1)
% hold on
% line([X(1) X(3)],[X(2) X(4)])


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


function sector = calculate_circle_sector(center, radius, point_a, point_b)
    center_to_a = point_a - center;
    center_to_b = point_b - center;

    cos_theta = dot(center_to_a, center_to_b) / (norm(center_to_a) * norm(center_to_b));
    theta= acos(cos_theta);
    % norm(point_a - center)
    % norm(point_b - center)
    % dot(norm(point_a - center), norm(point_b - center))
    % angle
    sector = theta * radius;
end

function rope_length = calculate_rope_length(circle_center_list,circle_radii_list,points_list)
    sum_of_circle_sector_length = 0;
    for i = 1:size(circle_center_list, 2)
        circle_sector_length = calculate_circle_sector(circle_center_list(:, i), circle_radii_list(i),points_list(:,2*i-1),points_list(:,2*i));
        sum_of_circle_sector_length = sum_of_circle_sector_length + circle_sector_length;
    end
    % circle_sector_length = calculate_circle_sector(circle_center_list(1), circle_radii_list(1),points_list(:,1),points_list(:,2))
    % circle_sector_length = calculate_circle_sector(circle_center_list(:, 2), circle_radii_list(2),points_list(:,3),points_list(:,4))   
    % circle_sector_length = calculate_circle_sector(circle_center_list(3), circle_radii_list(3),points_list(:,5),points_list(:,6))
    
    rope_segment_length = norm(points_list(:,1) - points_list(:,3));
    rope_segment_length = rope_segment_length + norm(points_list(:,4) - points_list(:,6));
    rope_segment_length = rope_segment_length + norm(points_list(:,2) - points_list(:,5));

    rope_length = sum_of_circle_sector_length + rope_segment_length;

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

        xList = [x];
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
        Xrot = x;
        iter = iterations;
        plot(xList(1:4:iterations * 4), ...
            xList(2:4:iterations * 4), '*', ...
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


function plot_line_from_matrix(mat)
    line(mat(1, :), mat(2, :), 'Color', 'red')
end

function drawSphere(center, radius, fmtStr)
    % https://www.youtube.com/watch?v=_NAzGCYowkc&ab_channel=ANALOGTecHII
    t = 0:pi/100:2 * pi;
    x = cos(t) * radius + center(1);
    y = sin(t) * radius + center(2);
    patch(x,y, 'blue')
% fplot(@(x) center(1) + radius * cos(x), @(y) center(2) + radius * sin(y), fmtStr)
end
