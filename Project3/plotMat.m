fidX = fopen('X.dat')
fidY = fopen('Y.dat')
rows = 300 * 365
X = fscanf(fidX,'%g',[10 rows]).'
Y = fscanf(fidY, '%g', [10 rows]).'
plot(X,Y)
%legend('Solar system Simulation');