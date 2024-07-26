function [G]=MakeAngles(firstCoordinate)
n=size(firstCoordinate,1);
a=ones(n,1)*[1:n];
a=a(:);
b=(ones(n,1)*[1:n])';
b=b(:);
c=(ones(n,1)*[2:n+1])';
c=c(:);
D=[b,a,c];
[r,c]=find(D==n+1);
D(r,:)=[];
D([1:n:n^2-n],:)=[];
D(n-1:n-1:end,3)=D(n-1:n-1:end,3)+1;
D(end,:)=[];
G=[ang(D(:,2),D(:,3),D(:,1),firstCoordinate(:,1),firstCoordinate(:,2))].*180/pi;
G=[D,G];
end


