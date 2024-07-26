clc;clear;chisquare_val =  2.4477;
format long g
CHOICE=menu('Select The Project You Want TO see','Design Geodetic Networks','FOD xp yp','FOD6','END');
while (CHOICE~=4)
    a=[];A=[];X=[];H1=[];
    switch CHOICE
        case 1
            NormQx_Ccs=[0;0;0;0];NormQx=NormQx_Ccs;rbar=NormQx_Ccs;chisquare_val =  2.4477;
            %----------------------------Add_Data-----------------------------------------------------------------
            [n1,p1]=uigetfile('*.txt','firstCoordinate.txt');
            firstCoordinate=strcat(p1,n1);
            firstCoordinate=dlmread(firstCoordinate);
            n=size(firstCoordinate,1);
            choiceee=0;
            while (choiceee~=5)
                choiceee=menu('Inpute Your Data','Lines','Angles','Azimuth','No Data','Continue');
                switch choiceee
                    case 1
                        [n1,p1]=uigetfile('*.txt','Lines.txt');
                        Lines=strcat(p1,n1);
                        Lines=dlmread(Lines);
                        errorline=input('errorline='); %0.0005
                    case 2
                        [n1,p1]=uigetfile('*.txt','Angles.txt');
                        Angles=strcat(p1,n1);
                        Angles=dlmread(Angles);
                        errorangle=input('errorangle='); %0.9s
                        errorangle=errorangle/3600;
                    case 3
                        [n1,p1]=uigetfile('*.txt','firstCoordinate.txt');
                        Azimuth=strcat(p1,n1);
                        Azimuth=dlmread(Azimuth);
                        errorAzimuth=input('errorAzimuth=');%1s
                        
                        errorAzimuth=errorAzimuth/3600;
                    case 4
                        KKK=1;
                        for i=1:n
                            
                            for j=i+1:n
                                LLL(KKK,:)=[i,j];
                                KKK=KKK+1;
                            end
                        end
                        errorline=input('errorline=');
                        d=(sum((firstCoordinate(LLL(:,1),:)-firstCoordinate(LLL(:,2),:)).^2,2)).^0.5;
                        d=d+(2*errorline*rand(KKK-1,1)-errorline.*ones(KKK-1,1));
                        Lines=[LLL,d];
                        errorangle=input('errorangle='); %1s
                        errorAzimuth=input('errorAzimuth=');%1s
                        errorangle=errorangle/3600;
                        errorAzimuth=errorAzimuth/3600;
                        Angles=MakeAngles(firstCoordinate);
                        Angles(:,4)= Angles(:,4)+(2*errorangle*rand(size(Angles,1),1)-errorangle.*ones(size(Angles,1),1));
                        Azimuth=MakeAzimuth(LLL,firstCoordinate(:,1),firstCoordinate(:,2));
                        Azimuth=Azimuth+(2*errorAzimuth*rand(size(Azimuth,1),1)-errorAzimuth.*ones(size(Azimuth,1),1));
                        Azimuth=[LLL,Azimuth];
                        dlmwrite('Azimuth.txt',Azimuth,'newline','pc','precision',12);
                        dlmwrite('Lines.txt',Lines,'newline','pc','precision',12);
                        dlmwrite('Angles.txt',Angles,'newline','pc','precision',12);
                end
            end
            %----------------------------------------------------------------------------------------------------
            
            %                         Lines=dlmread('C:\Users\ANDISHE\Desktop\Data\Lines.txt');
            %                         Angles=dlmread('C:\Users\ANDISHE\Desktop\Data\Angles.txt'); %Ftom,Center,To
            %                         firstCoordinate=dlmread('C:\Users\ANDISHE\Desktop\Data\firstCoordinate.txt');
            %                         errorangle=0.9/3600;
            %                         errorline=0.0005;
            %                         optimumellipsesize_Qc=0.005;
            %----------------------------------------------------------------------------------------------------
            optimumellipsesize_Qc=input('optimumellipsesize_Qc='); %0.005
            choice=menu('select Type of your data you want to use','Lines & Angles','Just Lines','Just Angles','Line & Azimuth','Continue');
            while (choice~=5)
                yy=[];Q=[];L=[];angles=[];yo=[];a=[];A=[];D=[];x0=[];delx=[];yzero=[];sxyRE=[];S=[];
                n=size(firstCoordinate,1);Qc=(optimumellipsesize_Qc/chisquare_val)^2.*eye(2*n);
                X = sym('x',[1 n]);Y = sym('y',[1 n]);
                if (choice==1) %Lines%angles
                    yy=[Lines(:,3);Angles(:,4)*pi/180];
                    Q=diag([errorline*ones(1,size(Lines,1)) errorangle*pi/(180)*ones(1,size(Angles,1))]).^2 ;
                    L=[dist(Lines(:,1),Lines(:,2),X',Y')];
                    angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)];
                    yo=[L;angles];
                else if (choice==2)  %%Lines
                        Q=diag([errorline*ones(1,size(Lines,1))]).^2 ;
                        yy=[Lines(:,3)];
                        L=[dist(Lines(:,1),Lines(:,2),X',Y')];
                        yo=[L];
                    else if (choice==3) %% Angles
                            Q=diag([errorangle*pi/(180)*ones(1,size(Angles,1))]).^2 ;
                            yy=[Angles(:,4)*pi/180];
                            angles=[gizmt(Angles(:,2),Angles(:,3),Angles(:,1),X,Y)]; %angles
                            yo=[angles];
                        else if (choice==4) %% Lines & Azimuth
                                Q=diag([errorline*ones(1,size(Lines,1)),errorAzimuth*pi/(180)*ones(1,size(Azimuth,1))]).^2 ;
                                yy=[Lines(:,3);Azimuth(:,3)*pi/180];
                                L=[dist(Lines(:,1),Lines(:,2),X',Y')];
                                azimuth=[AzImUtH(Azimuth(:,1:2),X',Y')];
                                yo=[L;azimuth];
                            end
                        end
                    end
                end
                W=inv(Q);
                non(1:2:2*n)=X;non(2:2:2*n)=Y;
                a=jacobian(yo,non);
                X_coordinate=firstCoordinate(:,1);
                Y_coordinate=firstCoordinate(:,2);
                x0=double(subs(non,[X Y],[X_coordinate' Y_coordinate']))';
                D(1,1:2:2*n)=mean(firstCoordinate(:,1))*ones(1,n);
                D(2,2:2:2*n)=mean(firstCoordinate(:,2))*ones(1,n);
                delx=1;repitition=0;
                while (norm(delx)> 1e-8)
                    A=double(subs(a,non,x0'));
                    
                    if (choice==1)
                        yzero=[double(subs(L,non,x0'));ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
                        
                        D(3,2:2:2*n)=X_coordinate';
                        D(3,1:2:2*n)=-Y_coordinate';
                    else if (choice==2)
                            yzero=[double(subs(L,non,x0'))];clc;
                            D(3,2:2:2*n)=X_coordinate';
                            D(3,1:2:2*n)=-Y_coordinate';
                            
                        else if (choice==3)
                                yzero=[ang(Angles(:,2),Angles(:,3),Angles(:,1),X_coordinate,Y_coordinate)];clc
                                
                                D(3,2:2:2*n)=X_coordinate';
                                D(3,1:2:2*n)=-Y_coordinate';
                                D(4,1:2:2*n)=X_coordinate';
                                D(4,2:2:2*n)=Y_coordinate';
                            else if (choice==4)
                                    yzero=[double(subs(L,non,x0'));MakeAzimuth([Azimuth(:,1),Azimuth(:,2)],X_coordinate,Y_coordinate)*pi/180];clc
                                end
                            end
                        end
                    end
                    delx=inv(A'*W*A+D'*D)*A'*W*(yy-yzero);
                    x0=x0+delx;
                    X_coordinate=x0(1:2:2*n);
                    Y_coordinate=x0(2:2:2*n);
                    repitition=repitition+1
                end;
                ehat=(yy-yzero);H=D;H1=H;
                yhat=yzero;
                Qxhat=inv(A'*W*A+D'*D)-D'*inv(D*D'*D*D')*D;
                df=size(yy,1)-rank(A);
                Q2=ehat'*W*ehat/df;rbar(choice,1)=df/size(yy,1);
                Qxhat=Qxhat.*Q2;
                f=diag(Qxhat);
                sxy=[f(1:2:2*n) f(2:2:2*n) diag(Qxhat(1:2:2*n,2:2:2*n))];
                el=[sqrt(0.5*(sxy(:,1)+sxy(:,2)+sqrt((sxy(:,1)-sxy(:,2)).^2+4.*sxy(:,3).^2))),sqrt(0.5*(sxy(:,1)+sxy(:,2)-sqrt((sxy(:,1)-sxy(:,2)).^2+4.*sxy(:,3).^2))),atan(2*sxy(:,3)./(sxy(:,1)-sxy(:,2)))];
                % sxyRE=cov{xi yi xj yj xiyi xjyj xixj yiyk xiyj xjyi XYcoordinate_i XYcoordinate_j i&j}
                k=1;
                for i=1:n
                    for j=i+1:n
                        sxyRE(k,:)=[sxy(i,1) , sxy(i,2) , sxy(j,1) , sxy(j,2) , sxy(i,3) , sxy(j,3) , Qxhat(2*i-1,2*j-1), Qxhat(2*i,2*j) , Qxhat(2*i-1,2*j) , Qxhat(2*j-1,2*i) , X_coordinate(i) , Y_coordinate(i) , X_coordinate(j) , Y_coordinate(j) , i , j];
                        k=k+1;
                    end
                end
                thetaa=atan(  (sxyRE(:,14)-sxyRE(:,12))./(sxyRE(:,13)-sxyRE(:,11)) );
                % {s2dxij s2dyij sdxdyij}
                sxyRE=[sxyRE(:,1)+sxyRE(:,3)-2*sxyRE(:,7) , sxyRE(:,2)+sxyRE(:,4)-2*sxyRE(:,8) , sxyRE(:,5)+sxyRE(:,6)-(sxyRE(:,9)+sxyRE(:,10)) , 0.5*(sxyRE(:,11)+sxyRE(:,13)) , 0.5*(sxyRE(:,12)+sxyRE(:,14)) , sxyRE(:,15) , sxyRE(:,16)];
                elRelative=[sqrt(0.5*(sxyRE(:,1)+sxyRE(:,2)+sqrt((sxyRE(:,1)-sxyRE(:,2)).^2+4.*sxyRE(:,3).^2))),sqrt(0.5*(sxyRE(:,1)+sxyRE(:,2)-sqrt((sxyRE(:,1)-sxyRE(:,2)).^2+4.*sxyRE(:,3).^2))),atan(2*sxyRE(:,3)./(sxyRE(:,1)-sxyRE(:,2)))];%a,b,theta
                %%plotEllipse
                if (choice==1)
                    scale=1;
                    plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Angles',[Lines(:,1),Lines(:,2)]);
                    ShowAngles(3,[X_coordinate,Y_coordinate],Angles(:,2),Angles(:,1),Angles(:,3));
                else if (choice==2)
                        scale=1;
                        plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines',[Lines(:,1),Lines(:,2)]);
                    else if (choice==3)
                            scale=1;
                            plotEllipseRotatedd(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Angles');
                        else if (choice==4)
                                scale=1;
                                plotEllipseRotated(chisquare_val*el(:,1),chisquare_val*el(:,2),[X_coordinate,Y_coordinate],wrapTo2Pi(el(:,3)),scale,'Lines & Azimuth',[Lines(:,1),Lines(:,2);Azimuth(:,1),Azimuth(:,2)]);
                            end
                        end
                    end
                end
                scale=1;
                plotEllipseRilative(chisquare_val*elRelative(:,1),chisquare_val*elRelative(:,2),[sxyRE(:,4:5)],wrapTo2Pi(elRelative(:,3)+thetaa),[sxyRE(:,6:7)],scale);
                S_Qcs=eye(2*n)-D'*inv(D*D')*D;
                Qcs=S_Qcs*Qc*Qc';
                NormQx(choice,1)=norm(Qxhat);
                NormQx_Ccs(choice,1)=norm(Qxhat-Qcs);
                Ri=diag(eye(size(yy,1))-A*inv(A'*W*A+H1'*H1)*A'*W);
                Ri_min_max(choice,:)=[find(Ri==min(Ri)),min(Ri),find(Ri==max(Ri)),max(Ri)]
                choice=menu('select Type of your data you want to use','Lines & Angles','Just Lines','Just Angles','Line & Azimuth','Continue');
            end
            char(string({'rbar_Lines&Angles= ';'rbar_Lines= ';'rbar_Angles= ';'rbar_Lines&Azimuth= '})+string(rbar))
            char(string({'normQx_Lines&Angles= ';'normQx_Lines= ';'normQx_Angles= ';'normQx_Lines&Azimuth= '})+string(NormQx))
            char(string({'NormQx_Ccs_Lines&Angles= ';'NormQx_Ccs_Lines= ';'NormQx_Ccs_Angles= ';'NormQx_Ccs_Lines&Azimuth= '})+string(NormQx_Ccs))
            choiceZoD=menu('ZOD To find the best point that Can be fixed','Continue');
            % ZOD
            H(3:end,:)=[];
            for i=1:n
                D=[];D=zeros(size(H));
                D(1,2*i-1)=20000;
                D(2,2*i)=12000;
                S_ZOD=eye(2*n)-H'*inv(D*H')*D;
                QXZOD=S_ZOD*Qxhat*S_ZOD';
                NormQx_ZOD(i,1)=norm(QXZOD);
            end
            i=find(NormQx_ZOD==min(NormQx_ZOD))
            choiceSOD=menu('FOD','Continue');
            % FOD
            X_coordinate_Adjustment=X_coordinate;
            Y_coordinate_Adjustment=Y_coordinate;
            X_Y_coordinate_Adjustment=[X_coordinate_Adjustment,Y_coordinate_Adjustment];
            min_ri(1)=min(Ri);
            ki=1;
            CYC=1;
            min_ri=min_ri(end);
            while(CYC~=2)
                POintNumber=input('POintNumber=');
                N=size(POintNumber,2);
                disp('input size of your box: ');
                m_L=input('x_Left=');
                if (N~=size(m_L,2))
                    m_L=m_L.*ones(1,N);
                end
                m_R=input('x_Right=');
                if (N~=size(m_R,2))
                    m_R=m_R.*ones(1,N);
                end
                n_L=input('y_Down=');
                if (N~=size(n_L,2))
                    n_L=n_L.*ones(1,N);
                end
                n_R=input('y_UP=');
                if (N~=size(n_R,2))
                    n_R=n_R.*ones(1,N);
                end
                precision=input('precision=');
                if (N~=size(precision,2))
                    precision=precision.*ones(1,N);
                end
                F_step=input('First_step=');
                if (N~=size(F_step,2))
                    F_step=F_step.*ones(1,N);
                end
                for PN=1:N
                    i=POintNumber(PN);pt=i;
                    XX=X_coordinate(i);
                    YY=Y_coordinate(i);
                    syms x y
                    A_FOD_i=subs(a,[X(i),Y(i)],[x,y]);
                    if (i==1)
                        X_FOD_i=[X(2:end)];Y_FOD_i=[Y(2:end)];
                        X_coordinate_FOD_i=[X_coordinate(2:end)'];Y_coordinate_FOD_i=[Y_coordinate(2:end)'];
                    else if (i==n)
                            X_FOD_i=[X(1:end-1)];Y_FOD_i=[Y(1:end-1)];
                            X_coordinate_FOD_i=[X_coordinate(1:end-1)'];Y_coordinate_FOD_i=[Y_coordinate(1:end-1)'];
                        else
                            X_FOD_i=[X(1:i-1),X(i+1:end)];Y_FOD_i=[Y(1:i-1),Y(i+1:end)];
                            X_coordinate_FOD_i=[X_coordinate(1:i-1)',X_coordinate(i+1:end)'];Y_coordinate_FOD_i=[Y_coordinate(1:i-1)',Y_coordinate(i+1:end)'];
                        end
                    end
                    A_FOD_i=subs(A_FOD_i,[X_FOD_i,Y_FOD_i],[X_coordinate_FOD_i Y_coordinate_FOD_i]);
                    ymin=YY-n_L(PN);ymax=YY+n_R(PN);
                    xmin=XX-m_L(PN);xmax=XX+m_R(PN);
                    step_X=(xmax-xmin)/F_step(PN);
                    step_Y=(ymax-ymin)/F_step(PN);
                    tic;
                    while (((step_X<precision(PN))&&(step_Y<precision(PN)))==0)
                        KK=0;
                        for x=xmin:step_X:xmax
                            K=1;KK=KK+1;
                            for y=ymin:step_Y:ymax
                                a_FOD_i=eval(A_FOD_i);
                                %RR(KK,K)=min(diag(eye(size(yy,1))-a_FOD_i*inv(a_FOD_i'*W*a_FOD_i+H1'*H1)*a_FOD_i'*W));K=K+1;
                                RR(KK,K)=min(diag(eye(size(yy,1))-a_FOD_i*inv(a_FOD_i'*a_FOD_i+H1'*H1)*a_FOD_i'));K=K+1;
                            end
                        end
                        [row,column]=find(RR==max(max(RR)));
                        row=row(1);
                        column=column(1);
                        ymax=ymin+(column)*step_Y;
                        if (ymax > YY+n_R(PN))
                            ymax=YY+n_R(PN);
                        end
                        ymin=ymin+(column-2)*step_Y;
                        if (ymin < YY-n_L(PN))
                            ymin=YY-n_L(PN);
                        end
                        xmax=xmin+(row)*step_X;
                        if (xmax > XX+m_R(PN))
                            xmax=XX+m_R(PN);
                        end
                        xmin=xmin+(row-2)*step_X;
                        if (xmin < XX-m_L(PN))
                            xmin=XX-m_L(PN);
                        end
                        step_X=2*step_X/F_step(PN);
                        step_Y=2*step_Y/F_step(PN);
                        R=RR;
                        RR=[];
                    end
                    toc
                    Yi=ymin+(column-1)*step_Y;
                    Xi=xmin+(row-1)*step_X;
                    ki=ki+1;
                    min_ri(ki)=max(max(R));
                    if  (min_ri(ki)> min_ri(ki-1))
                        X_coordinate(pt)=Xi;
                        Y_coordinate(pt)=Yi;
                    end
                    
                end
                disp('min_ri')
                disp(min_ri)
                disp('min_ri_Last')
                disp(min_ri(ki))
                CYC=menu('FOD','Continue','end');
            end
            
            
            
            
            hold on
            axis equal;
            for i=1:size(X_coordinate,1)
                plot(X_coordinate(i,1),Y_coordinate(i,1),'.');
                r=char(string(' P')+string(i));
                text(X_coordinate(i,1),Y_coordinate(i,1),r,'Color','blue','FontSize',8);
                for j=i+1:size(X_coordinate,1);
                    line([X_coordinate(i,1);X_coordinate(j,1)],[Y_coordinate(i,1);Y_coordinate(j,1)]);
                end
            end
            X_Y_cs_after_FOD=[X_coordinate,Y_coordinate];
            disp('X_Y_cs_after_FOD - X_Y_coordinate_Adjustment')
            disp(X_Y_cs_after_FOD-X_Y_coordinate_Adjustment)
            dlmwrite('X_Y_cs_after_FOD.txt',X_Y_cs_after_FOD,'newline','pc','precision',12);
            
            
            
            
        case 2
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
                    y1=-5.2255;y2=99;n=0.5;
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
        case 3
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
            plot(X,Y,'*r')
    end
    CHOICE=menu('Select The Project You Want TO see','Design Geodetic Networks','FOD xp yp','FOD6','END');
end