fidX = fopen('X.dat')
fidY = fopen('Y.dat')
X = fscanf(fidX,'%g',[10 1000]).'
Y = fscanf(fidY, '%g', [10 1000]).'
plot(X,Y)
%legend('Solar system Simulation');