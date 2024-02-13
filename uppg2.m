a = [-1; 1.5];
b = [3; 0.5];
c = [0; -2];

ra = 1;
rb = 1.2;
rc = 1.7;

[l, X] = langd(ra, rb, rc, a, b, c);

disp("Längd");
disp(l);

axis equal
hold on
for i = [1:3]
    line(X(1:2:end,i), X(2:2:end,i));
end

drawSphere(a, ra);
drawSphere(b, rb);
drawSphere(c, rc);

%% Experimentell störningsanalys


% Kod för uppgift 2...
function [L,X]=langd(ra,rb,rc,a,b,c)
    a_b_guess = get_start_guesses(a, b, c);
    a_b_coords= punkter(a_b_guess, ra, rb, a, b);
    
    a_c_guess = get_start_guesses(a, c, b);
    a_c_coords = punkter(a_c_guess, ra, rc, a, c);


    b_c_guess = get_start_guesses(b, c, a);
    b_c_coords = punkter(b_c_guess, rb, rc, b, c);

    points_list = [a_b_coords(1:2), a_c_coords(1:2), a_b_coords(3:4), b_c_coords(1:2), a_c_coords(3:4), b_c_coords(3:4)];
    L = calculate_rope_length([a, b, c], [ra, rb, rc], points_list);

    X =  [a_b_coords, b_c_coords, a_c_coords];
end


function guesses = get_start_guesses(from, to, other)
    from_to = from + (from - other);
    to_from = to + (to - other);

    guesses = [from_to; to_from];
end


function sector = calculate_circle_sector(center, radius, point_a, point_b)
    center_to_a = point_a - center;
    center_to_b = point_b - center;

    cos_theta = dot(center_to_a, center_to_b) / (norm(center_to_a) * norm(center_to_b));
    theta= acos(cos_theta);
    sector = theta * radius;
end

function rope_length = calculate_rope_length(circle_center_list,circle_radii_list,points_list)
    sum_of_circle_sector_length = 0;
    for i = 1:size(circle_center_list, 2)
        circle_sector_length = calculate_circle_sector(circle_center_list(:, i), circle_radii_list(i),points_list(:,2*i-1),points_list(:,2*i));
        sum_of_circle_sector_length = sum_of_circle_sector_length + circle_sector_length;
    end
    
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

    F = @(x) [(x(1)-a1)^2+(x(2)-a2)^2-ra^2;
              (x(3)-b1)^2+(x(4)-b2)^2-rb^2;
              (x(1)-x(3))*(x(1)-a1)+(x(2)-x(4))*(x(2)-a2);
              (x(1)-x(3))*(x(3)-b1)+(x(2)-x(4))*(x(4)-b2)];

    jacobian = @(x) [2*(x(1)-a1) 2*(x(2)-a2) 0 0;
                     0 0 2*(x(3)-b1) 2*(x(4) -b2);
                     2*x(1)-a1-x(3) 2*x(2)-a2-x(4) a1-x(1) a2-x(2);
                     x(3)-b1 x(4)-b2 x(1)+b1-2*x(3) x(2)+b2-2*x(4)];

    x = X0;

    tolerance = 1e-10;
    iterations = 0;

    xPrev = [inf inf inf inf]'; % x + 2 * tolerance;

    while max(abs(x - xPrev)) > tolerance
        xPrev = x;
        x = x - jacobian(x) \ F(x);
        iterations = iterations + 1;
    end

    Xrot = x;
    iter = iterations;
end

%% Drawing functions
function drawSphere(center, radius)
    % https://www.youtube.com/watch?v=_NAzGCYowkc&ab_channel=ANALOGTecHII
    t = 0:pi/100:2 * pi;
    x = cos(t) * radius + center(1);
    y = sin(t) * radius + center(2);
    patch(x,y, 'blue')
end
