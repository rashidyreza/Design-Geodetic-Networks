clc;clear;%char(string('normQ= ')+string(A))
Lines=dlmread('C:\Users\ANDISHE\Desktop\Data\Lines.txt');
Angles=dlmread('C:\Users\ANDISHE\Desktop\Data\Angles.txt'); %Ftom,Center,To
firstCoordinate=dlmread('C:\Users\ANDISHE\Desktop\Data\firstCoordinate.txt');




choice=0;
    choice=menu('select Type of your data you want to use','Lines & Angles','Just Lines','Just Angles');
yy=[];Q=[];L=[];angles=[];yo=[];a=[];A=[];D=[];x0=[];delx=[];yzero=[];sxyRE=[];




n=size(firstCoordinate,1);
X = sym('x',[1 n]);Y = sym('y',[1 n]);
if (choice==1) %Lines%angles
    yy=[Lines(:,3);Angles(:,4)*pi/180];
    Q=diag([0.005*ones(1,size(Lines,1)) 0.00025*pi/(180)*ones(1,size(Angles,1))]).^2 ;
    L=[dist(Lines(:,1),Lines(:,2),X',Y')];
    angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)];
    yo=[L;angles];
else if (choice==2)  %%Lines
        Q=diag([0.005*ones(1,size(Lines,1))]).^2 ;
        yy=[Lines(:,3)];
        L=[dist(Lines(:,1),Lines(:,2),X',Y')];
        yo=[L];
    else if (choice==3) %% Angles
            Q=diag([0.00025*pi/(180)*ones(1,size(Angles,1))]).^2 ;
            yy=[Angles(:,4)*pi/180];
            angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)]; %angles
            yo=[angles];
        end
    end
end
W=inv(Q);
non(1:2:2*n)=X;non(2:2:2*n)=Y;
a=jacobian(yo,non);
X_coordinate=firstCoordinate(:,1);
Y_coordinate=firstCoordinate(:,2);
x0=double(subs(non,[X Y],[X_coordinate' Y_coordinate']))';
D(1,1:2:2*n)=20000*ones(1,n);
D(2,2:2:2*n)=12000*ones(1,n);
delx=1;repitition=0;
while (norm(delx)> 1e-8)
A=double(subs(a,non,x0'));
D(3,2:2:2*n)=X_coordinate';
D(3,1:2:2*n)=-Y_coordinate';
if (choice==1)
    yzero=[double(subs(L,non,x0'));ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
else if (choice==2)
        yzero=[double(subs(L,non,x0'))];clc;
    else if (choice==3)
            yzero=[ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
            D(4,1:2:2*n)=X_coordinate';
            D(4,2:2:2*n)=Y_coordinate';
        end
    end
end
delx=inv(A'*W*A+D'*D)*A'*W*(yy-yzero);
x0=x0+delx;
X_coordinate=x0(1:2:2*n);
Y_coordinate=x0(2:2:2*n);
repitition=repitition+1
end;
ehat=(yy-yzero);
yhat=yzero;
Qxhat=inv(A'*W*A+D'*D)-D'*inv(D*D'*D*D')*D;
df=size(yy,1)-rank(A);
Q2=ehat'*W*ehat/df;
Qxhat=Qxhat.*Q2;QxAngles=trace(Qxhat);
f=diag(Qxhat);
sxy=[f(1:2:23) f(2:2:24) diag(Qxhat(1:2:23,2:2:24))];
el=[sqrt(0.5*(sxy(:,1)+sxy(:,2)+sqrt((sxy(:,1)-sxy(:,2)).^2+4.*sxy(:,3).^2))),sqrt(0.5*(sxy(:,1)+sxy(:,2)-sqrt((sxy(:,1)-sxy(:,2)).^2+4.*sxy(:,3).^2))),atan(2*sxy(:,3)./(sxy(:,1)-sxy(:,2)))];
% sxyRE=cov{xi yi xj yj xiyi xjyj xixj yiyk xiyj xjyi XYcoordinate_i XYcoordinate_j i&j}
k=1;
for i=1:n
    for j=i+1:n
     sxyRE(k,:)=[sxy(i,1) , sxy(i,2) , sxy(j,1) , sxy(j,2) , sxy(i,3) , sxy(j,3) , Qxhat(2*i-1,2*j-1), Qxhat(2*i,2*j) , Qxhat(2*i-1,2*j) , Qxhat(2*j-1,2*i) , X_coordinate(i) , Y_coordinate(i) , X_coordinate(j) , Y_coordinate(j) , i , j];
     k=k+1;
    end
end
% {s2dxij s2dyij sdxdyij}
sxyRE=[sxyRE(:,1)+sxyRE(:,3)-2*sxyRE(:,7) , sxyRE(:,2)+sxyRE(:,4)-2*sxyRE(:,8) , sxyRE(:,5)+sxyRE(:,6)-(sxyRE(:,9)+sxyRE(:,10)) , 0.5*(sxyRE(:,11)+sxyRE(:,13)) , 0.5*(sxyRE(:,12)+sxyRE(:,14)) , sxyRE(:,15) , sxyRE(:,16)]  
elRelative=[sqrt(0.5*(sxyRE(:,1)+sxyRE(:,2)+sqrt((sxyRE(:,1)-sxyRE(:,2)).^2+4.*sxyRE(:,3).^2))),sqrt(0.5*(sxyRE(:,1)+sxyRE(:,2)-sqrt((sxyRE(:,1)-sxyRE(:,2)).^2+4.*sxyRE(:,3).^2))),atan(2*sxyRE(:,3)./(sxyRE(:,1)-sxyRE(:,2)))];%a,b,theta

chisquare_val =  2.4477;%%plotEllipse
if (choice==1)
    scale=20000;
    plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Angles',[Lines(:,1),Lines(:,2)]);
    ShowAngles(3,[X_coordinate,Y_coordinate],Angles(:,2),Angles(:,1),Angles(:,3));
else if (choice==2)
        scale=2000;
        plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines',[Lines(:,1),Lines(:,2)]);
    else if (choice==3)
            scale=20000;
            plotEllipseRotatedd(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Angles');
        end
    end
end
scale=2000;
plotEllipseRilative(chisquare_val*elRelative(:,1),chisquare_val*elRelative(:,2),[sxyRE(:,4:5)],wrapTo2Pi(elRelative(:,3)),[sxyRE(:,6:7)],scale);

