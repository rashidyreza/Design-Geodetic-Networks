function ShowAngles(R,C,Centre,from,To)
hold on;R1=R;
N = 100;
for i=1:size(Centre,1);
    theta2 =gizm(C(Centre(i),1),C(Centre(i),2),C(from(i),1),C(from(i),2));
    theta1 =gizm(C(Centre(i),1),C(Centre(i),2),C(To(i),1),C(To(i),2));
    theta=linspace(theta1,theta2,N);
    X = R1*sin( theta );
    Y = R1*cos( theta );
    X = X + C(Centre(i),1)*ones(1,N);
    Y = Y + C(Centre(i),2).*ones(1,N);
    % Plot     %----------------------------------------
    plot(X,Y);R1=R+6*rand(1,1);
end
end
