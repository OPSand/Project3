fidX = fopen('X.dat')
fidY = fopen('Y.dat')
rows = 300 * 365
X = fscanf(fidX,'%g',[10 rows]).';
Y = fscanf(fidY, '%g', [10 rows]).';
%plot(X,Y)
for i=1:10
   plot(X(:,i),Y(:,i),'color',rand(1,3))
   hold on
end
legend('Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptun','Pluto',10)
