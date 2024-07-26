function [G]=gizmt(Centre,from,To,X,Y)
for i=1:size(To,1)
G(i,1)=atan((X(Centre(i))-X(To(i)))/(Y(Centre(i))-Y(To(i))))-atan((X(Centre(i))-X(from(i)))/(Y(Centre(i))-Y(from(i))));
end
