function ObjVal = Schaffer(Chrom,switch1);
    %Dim=size(Chrom,2);%返回矩阵的列数
    %[Nind, Nvar] = size(Chrom);
    [Row,Dim] = size(Chrom);
    value = zeros(1,Row);
   for i=1:Row
    for j=1:Dim
        value(i) =value(i)+ Chrom(i,j).^2;
    end;
   end;
    ObjVal = 0.5 + (sin(sqrt(value)).^2 - 0.5)./(1+0.001.*value).^2;
