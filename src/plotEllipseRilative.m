function plotEllipseRilative(a,b,C,angle,KK,scale)
a=a*scale;b=b*scale;
axis equal;
% range to plot over
%------------------------------------
N = 100;hold on;
theta = 0:1/N:2*pi+1/N;
% Parametric equation of the ellipse
kkk=char(string('Relative Ellipse scale is ')+string(scale));
ylabel(kkk)
for i=1:size(angle,1)
    phi = angle(i);
    
    %----------------------------------------
    ellipse_x_r  = a(i)*cos( theta );
    ellipse_y_r  = b(i)*sin( theta );
    
    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    %let's rotate the ellipse to some angle phi
    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
    
    % Coordinate transform (since your ellipse is axis aligned)
    %----------------------------------------
    X = r_ellipse';
    X(1,:) = X(1,:) + C(i,1);
    X(2,:) = X(2,:) + C(i,2);
    % Plot     %----------------------------------------
    plot(X(1,:),X(2,:));
    r=char(string('P')+string(KK(i,1))+string(' P')+string(KK(i,2)));
    text(C(i,1),C(i,2),r,'Color','blue','FontSize',4);
     plot(C(i,1),C(i,2),'.g');

end

