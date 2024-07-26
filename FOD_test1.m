clc;clear;
choiceFOD=menu('FOD','x=50 y=?','x=? & y=?');
tic;
switch choiceFOD
    case 1
% y=50/sqrt(3)
syms x y 
A=[sqrt((x-0)^2+(y-0)^2);
    sqrt((x-100)^2+(y-0)^2);
    sqrt((x-50)^2+(y-100)^2)];
A=jacobian(A,[x,y]);x=50;
A=eval(A);
y1=0;y2=99;n=0.5;
for j=1:10
    K=1;
for y=y1:n:y2
    a=eval(A);
R(K)=min(diag(eye(3)-a*inv(a'*a)*a'));K=K+1;
end
i=find(R==max(R));
yy=[y1:n:y2];
Y=y1+(i-1)*n;
y1=y1+(i-2)*n;
y2=y1+(i)*n;
n=0.1*n;RR=R;
R=[];
figure;
plot(yy,RR)
xlabel('yp')
ylabel('min(ri)')
end
disp(Y)
Y-50/sqrt(3)

    case 2   
      syms x y 
A=[sqrt((x-0)^2+(y-0)^2);
    sqrt((x-100)^2+(y-0)^2);
    sqrt((x-50)^2+(y-100)^2)];
A=jacobian(A,[x,y]);
y1=25;y2=30;n=1;
x1=45;x2=55;
for j=1:5
    KK=0;
    for x=x1:n:x2
             K=1;KK=KK+1;
for y=y1:n:y2
    a=eval(A);
R(KK,K)=min(diag(eye(3)-a*inv(a'*a)*a'));K=K+1;
end
    end
    [row,column]=find(R==max(max(R)));
Y=y1+(column-1)*n;
X=x1+(row-1)*n;
y1=y1+(column-2)*n;
y2=y1+(column)*n;
x1=x1+(row-2)*n;
x2=x1+(row)*n;
n=0.1*n;
% RR=R;
R=[];
end
disp([Y,X])
Y-50/sqrt(3)  
end
toc