clc;clear;
Lines=dlmread('C:\Users\ANDISHE\Desktop\Data\Lines.txt');
Angles=dlmread('C:\Users\ANDISHE\Desktop\Data\Angles.txt'); %Ftom,Center,To
firstCoordinate=dlmread('C:\Users\ANDISHE\Desktop\Data\firstCoordinate.txt');
Q=diag([0.005*ones(1,size(Lines,1)) 0.00025*pi/(180)*ones(1,size(Angles,1))]).^2 ;
W=inv(Q);
yy=[Lines(:,3);Angles(:,4)*pi/180];
n=size(firstCoordinate,1);
X = sym('x',[1 n]);Y = sym('y',[1 n]);
L=[dist(Lines(:,1),Lines(:,2),X',Y')]; %Lines
angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)]; %angles
yo=[L;angles];
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
yzero=[double(subs(L,non,x0'));ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
D(3,2:2:2*n)=X_coordinate';
D(3,1:2:2*n)=-Y_coordinate';
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
chisquare_val =  2.4477;
scale=20000;
     plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Angles',[Lines(:,1),Lines(:,2)]);
     ShowAngles(1,[X_coordinate,Y_coordinate],Angles(:,2),Angles(:,1),Angles(:,3));
     %% Lines
Lines=dlmread('C:\Users\ANDISHE\Desktop\Data\Lines.txt');
firstCoordinate=dlmread('C:\Users\ANDISHE\Desktop\Data\firstCoordinate.txt');
Q=diag([0.005*ones(1,size(Lines,1))]).^2 ;
W=inv(Q);
yy=[Lines(:,3)];
n=size(firstCoordinate,1);
X = sym('x',[1 n]);Y = sym('y',[1 n]);
L=[dist(Lines(:,1),Lines(:,2),X',Y')]; %Lines
yo=[L];
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
yzero=[double(subs(L,non,x0'))];clc;
D(3,2:2:2*n)=X_coordinate';
D(3,1:2:2*n)=-Y_coordinate';
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
chisquare_val =  2.4477;
scale=2000;
     plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines',[Lines(:,1),Lines(:,2)]);
   
     %% Angles
Angles=dlmread('C:\Users\ANDISHE\Desktop\Data\Angles.txt'); %Ftom,Center,To
firstCoordinate=dlmread('C:\Users\ANDISHE\Desktop\Data\firstCoordinate.txt');
Q=diag([0.00025*pi/(180)*ones(1,size(Angles,1))]).^2 ;
W=inv(Q);
yy=[Angles(:,4)*pi/180];
n=size(firstCoordinate,1);
X = sym('x',[1 n]);Y = sym('y',[1 n]);
angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)]; %angles
yo=[angles];
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
yzero=[ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
D(3,2:2:2*n)=X_coordinate';
D(3,1:2:2*n)=-Y_coordinate';
D(4,1:2:2*n)=X_coordinate';
D(4,2:2:2*n)=Y_coordinate';
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
chisquare_val =  2.4477;
scale=20000;
plotEllipseRotatedd(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Angles');