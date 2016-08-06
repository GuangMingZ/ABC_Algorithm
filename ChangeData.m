clear;
clc;

Matri = load('E:\Schaffer\60-2result.txt');
[a,b] = size(Matri);
 
 S = transpose(Matri);
 S = reshape(S,a*b,1);
% for j=1:4
%     d = (j-1)*1000;
%  for i=800:870
%     S(i+d) = 0.00050;
%  end;
%  for i=870:900
%      S(i+d) = 0.00000;
%  end;
%   for i=900:1000
%      S(i+d) = 0.00000;
%   end;
% end;
 dlmwrite('Data.txt',S,'-append','delimiter', ' ');
 %save Data.txt -ascii S;