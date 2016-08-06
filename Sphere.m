
function ObjVal=Sphere(Foods)

%global Foods;
[a,b] = size(Foods);
value = zeros(1,a);
for i=1:a
    for j=1:b
    value(i) =value(i)+ Foods(i,j).^2;
    end;
end;
%S=Food.*Food;
ObjVal=value;%sum(S);

