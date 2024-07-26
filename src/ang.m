function [G]=ang(Centre,from,To,X,Y)
for i=1:size(To,1);
x1=X(Centre(i));
y1=Y(Centre(i));
x2=X(from(i));
y2=Y(from(i));
x3=X(To(i));
y3=Y(To(i));
G(i,1)=wrapTo2Pi(gizm(x1,y1,x3,y3)-gizm(x1,y1,x2,y2));
end