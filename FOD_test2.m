n=input('number of points= ');
R=500;
 k=input('Inpute Box Size= ');
theta=linspace(0,2*pi,n+1);theta=theta+pi/2;
X=R*cos(theta);Y=R*sin(theta);plot(X,Y);axis equal;
X=[X(1:n)';0];Y=[Y(1:n)';0];
X_c = sym('x',[1 n+1]);Y_c = sym('y',[1 n+1]);
for i=1:n+1
    r=char(string(' P')+string(i));
    text(X(i),Y(i),r,'Color','blue','FontSize',8);
end
KKK=1;
for i=1:n+1
  
    for j=i+1:n+1
        LLL(KKK,:)=[i,j];
        KKK=KKK+1;
    end
end
L=[dist(LLL(:,1),LLL(:,2),X_c,Y_c)];
A=jacobian(L,[X_c,Y_c]);

i=input('Input Number Of The Point you Want To Move: ');
XX=X(i);
YY=Y(i);
syms x y
A=subs(A,[X_c(i),Y_c(i)],[x,y]);
if (i==1)
    X_c=[X_c(2:end)];Y_c=[Y_c(2:end)];
    X=[X(2:end)'];Y=[Y(2:end)'];
else if (i==n+1)
        X_c=[X_c(1:end-1)];Y_c=[Y_c(1:end-1)];
        X=[X(1:end-1)'];Y=[Y(1:end-1)'];
    else
X_c=[X_c(1:i-1),X_c(i+1:end)];Y_c=[Y_c(1:i-1),Y_c(i+1:end)];
X=[X(1:i-1)',X(i+1:end)'];Y=[Y(1:i-1)',Y(i+1:end)'];

    end
end
A=subs(A,[X_c,Y_c],[X,Y]);
ymin=YY-k/2;ymax=YY+k/2;
xmin=XX-k/2;xmax=XX+k/2;
KK=0;h1=null(A);h1=h1';N=0.5*(n^2+n);
    for x=xmin:0.5:xmax
             K=1;KK=KK+1;
for y=ymin:0.5:ymax
    a=eval(A);H1=eval(h1);
R(KK,K)=min(diag(eye(N)-a*inv(a'*a+H1'*H1)*a'));K=K+1;
end
    end
    [row,column]=find(R==max(max(R)));
Y=ymin+(column-1)
X=xmin+(row-1)
hold on
plot(X,Y,'*')


