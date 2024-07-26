function [k]=dist(firstpoint,secondpoint,X,Y)
for i=1:size(firstpoint,1)
dx=X(secondpoint(i))-X(firstpoint(i));
dy=Y(secondpoint(i))-Y(firstpoint(i));
k(i,1)=sqrt(power(dx,2)+power(dy,2));
end

