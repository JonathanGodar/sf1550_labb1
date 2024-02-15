% Kod för uppgift 3...
%% a
load eiffel4.mat

trussplot(xnod,ynod,bars);
hold on
index = 1;
b = randn(size(A, 2),1);
b(2*index-1) = 1000;
tic
x = A\b;
toc

xbel = xnod + x(1:2:end); 
ybel = ynod + x(2:2:end);
disp([xbel(index),xnod(index)])
disp([ybel(index),ybel(index)])

plot(xnod(index), ynod(index), 'ro', 'MarkerSize', 3)
plot(xbel(index), ybel(index), 'bo', 'MarkerSize', 3)


trussplot(xbel,ybel,bars)

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
