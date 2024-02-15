% Kod för uppgift 3...
%% a
load eiffel4.mat

trussplot(xnod,ynod,bars);
hold on
index = 1;
b = zeros(size(A, 2),1);
b(2*index-1) = 1000;

x = A\b;


xbel = xnod + x(1:2:end); 
ybel = ynod + x(2:2:end);
disp([xbel(index),xnod(index)])
disp([ybel(index),ybel(index)])

plot(xnod(index), ynod(index), 'ro', 'MarkerSize', 3)
plot(xbel(index), ybel(index), 'bo', 'MarkerSize', 3)


trussplot(xbel,ybel,bars)
%%
load eiffel1.mat
e1Size = size(A,2);
e1Time = calculationTime(A,500);

load eiffel2.mat
e2Size = size(A,2);
e2Time = calculationTime(A,50);

load eiffel3.mat
e3Size = size(A,2);
e3Time = calculationTime(A,25);

load eiffel4.mat
e4Size = size(A,2);
e4Time = calculationTime(A,10);

sizeList = [e1Size, e2Size, e3Size, e4Size];
timeList = [e1Time, e2Time, e3Time, e4Time];
%%
loglog(sizeList,timeList)
grid on
hold on
plot(sizeList(1), timeList(1), 'ro', 'MarkerSize', 3)
plot(sizeList(2), timeList(2), 'ro', 'MarkerSize', 3)
plot(sizeList(3), timeList(3), 'ro', 'MarkerSize', 3)
plot(sizeList(4), timeList(4), 'ro', 'MarkerSize', 3)
x = linspace(sizeList(1), sizeList(4));
y = x - x(1) + timeList(1);
% plot(x,y)
% loglog(x,x)



%% Tidtabell

% Skapa en 4x4-matris T som innehåller beräkningstiderna.
% Raderna ska motsvara de olika modellerna (eiffel1-eiffel4) och
% kolumnerna de olika metoderna, ordnade som "Naiv", "LU",
% "Gles" och "Gles LU".

% Följande kod skapar en snygg tabell med resultaten:


tab=array2table(T,'VariableNames',{'Naiv' 'LU' 'Gles' 'Gles LU'},'RowNames',{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
disp(tab);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jmin,jmax]=kanslighet(A,metod)

%  Indata:
%
%  A     - matrisen
%  metod - villken metod som används:
%          1 = Naiv metod
%          2 = LU-faktorisering
%
%  Utdata:
%
%  jmin - index för minst känsliga nod
%  jmax - index för mest känsliga nod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

end

function time = calculationTime(A,iterations)
    b = randn(size(A, 2),1);
    tic
    for i = [1:iterations]
        x = A\b;
    end
    time = toc/iterations;
    
end
