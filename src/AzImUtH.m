function [G]=AzImUtH(LLL,X,Y)
Centre=LLL(:,1);
To=LLL(:,2);
for i=1:size(To,1)
    G(i,1)=atan((X(To(i))-X(Centre(i))) /(Y(To(i))-Y(Centre(i))));
end
end
