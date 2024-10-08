function [T] = FTCS(T,nx,ny,x1,x2,y1,y2,Fo,Bi,Ta)




    for i=1:nx
        for j=1:ny
           
            
            if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                continue
          
            elseif i==1 && j==1
                T(i,j)=2*Fo*(T(i+1,j)+T(i,j+1)+2*Bi*Ta)+(1-4*Fo-4*Bi*Fo)*T(i,j);
            elseif i==1 && j>1 && j<ny
                T(i,j)=Fo*(2*T(i+1,j)+T(i,j+1)+T(i,j-1)+2*Bi*Ta)+(1-4*Fo-2*Bi*Fo)*T(i,j);
                
            elseif j==ny && i==1
                T(i,j)=2*Fo*(T(i+1,j)+T(i,j-1)+2*Bi*Ta)+(1-4*Fo-4*Bi*Fo)*T(i,j);
            elseif i==nx && j==1
                T(i,j)=2*Fo*(T(i-1,j)+T(i,j+1)+2*Bi*Ta)+(1-4*Fo-4*Bi*Fo)*T(i,j);
            elseif i==nx && j==ny
                T(i,j)=2*Fo*(T(i-1,j)+T(i,j-1)+2*Bi*Ta)+(1-4*Fo-4*Bi*Fo)*T(i,j);
                
                
            elseif j==1 && i>1 && i<nx
                T(i,j)=Fo*(2*T(i,j+1)+T(i+1,j)+T(i-1,j)+2*Bi*Ta)+(1-4*Fo-2*Bi*Fo)*T(i,j);
            elseif j==ny && i>1 && i<nx
                T(i,j)=Fo*(2*T(i,j-1)+T(i+1,j)+T(i-1,j)+2*Bi*Ta)+(1-4*Fo-2*Bi*Fo)*T(i,j);
            
            elseif i==nx && j>1 && j<ny
                T(i,j)=Fo*(2*T(i-1,j)+T(i,j+1)+T(i,j-1)+2*Bi*Ta)+(1-4*Fo-2*Bi*Fo)*T(i,j);
            else
                T(i,j)=Fo*(T(i+1,j)+T(i-1,j)+T(i,j-1)+T(i,j+1)-4*T(i,j))+T(i,j);
            end
        end
    end
   

end

