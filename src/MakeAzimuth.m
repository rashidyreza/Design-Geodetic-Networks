function [G]=MakeAzimuth(LLL,X,Y)
Centre=LLL(:,1);
To=LLL(:,2);
for i=1:size(To,1)
    G(i,1)=wrapTo360(gizm(X(Centre(i)),Y(Centre(i)),X(To(i)),Y(To(i)))*180/pi);
end
end
