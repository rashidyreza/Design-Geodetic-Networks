function [G]=gizm(x1,y1,x2,y2)
dx=x2-x1;
dy=y2-y1;
v=atan(abs(x1-x2)/abs(y1-y2));
if (dx>0 && dy>0)
    G=v;
else if (dx>0 && dy<0)
        G=pi-v;
    else if (dx<0 && dy<0)
            G=pi+v;
        else if (dx<0 && dy>0)
                G=2*pi-v;
        
        else if(dx<0 && dy==0)
                G=3*pi/2;
            else if (dx>0 && dy==0)
                    G=pi/2;
                    else if (dx==0 && dy<0)
                                 G=pi;
                                 else if (dx==0 && dy>0)
                                         G=0;
                                     end
                        end

                end
            end
        end
    end
    end 
end