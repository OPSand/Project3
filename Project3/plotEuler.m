fidXeuler = fopen('Xeuler.dat')
fidYeuler = fopen('Yeuler.dat')
rows = 300 * 365
Xeuler = fscanf(fidXeuler,'%g',[10 rows]).';
Yeuler = fscanf(fidYeuler, '%g', [10 rows]).';
%plot(X,Y)
for i=1:10
   plot(Xeuler(:,i),Yeuler(:,i),'color',rand(1,3))
   hold on
end
legend('Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptun','Pluto',10)