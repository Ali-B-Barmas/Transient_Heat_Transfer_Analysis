function [Tm] = Lassonnen(T,nx,ny,x1,x2,y1,y2,Ta,Tin,w,Fo,Bi)





in=(2*(w+1)-1);
io=nx-in;
L=nx*ny-in*in;
CM=zeros(L);                             %Coefficient Matrix
GM=zeros(L,1);                           %Known Values Matrix
Tm=T;




A=1+4*Fo;
B=1+2*Fo*(2+Bi);
C=1+4*Fo*(1+Bi);
D=-Fo;

f=zeros(nx,ny);          
% f is a matrix that has number of each array witch will be used to create
% Coefficient Matrix
for j=1:ny
    for i=1:nx
        if (i>x1 && i<x2) && (j>y1 && j<y2)
                f(i,j)=nan;
        elseif j<y1
            f(i,j)=(j-1)*nx+i;
        elseif j>=y1 && j<=y2
            if i<x1
                f(i,j)=w*nx+(j-w-1)*io+( i );
            elseif i>x2
                f(i,j)=w*nx+(j-w-1)*io+( i-in );
            end
        elseif j>y2
            f(i,j)=w*nx+io*in+(j-(w+in)-1)*nx+i;
        end
    end
end         


    
    for j=1:y1-1
        for i=1:nx
        
         if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                continue
         
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         elseif j==1       
             
            if i==1
                GM(f(i,j))=Tm(i,j)+4*Bi*Fo*Ta;
                CM(f(i,j),f(i,j))=C;
                CM(f(i,j),f(i,j+1))=2*D;
                CM(f(i,j),f(i+1,j))=2*D;
            
            
            elseif i==nx
                GM(f(i,j))=Tm(i,j)+4*Bi*Fo*Ta;
                CM(f(i,j),f(i,j))=C;
                CM(f(i,j),f(i-1,j))=2*D;
                CM(f(i,j),f(i,j+1))=2*D;
            else 
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j+1))=2*D;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
         elseif j==y1-1
             if i==1
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;
                CM(f(i,j),f(i+1,j))=2*D;
                
             elseif i==nx
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;
                CM(f(i,j),f(i-1,j))=2*D;
                
             elseif i>=x1 && i<=x2
                GM(f(i,j))=Tm(i,j)-D*Tin;
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i,j-1))=D;
             else
                GM(f(i,j))=Tm(i,j);       
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;
         
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         else
             if i==1
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i+1,j))=2*D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;        
             elseif i==nx
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i-1,j))=2*D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;
             else
                GM(f(i,j))=Tm(i,j);       
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i,j-1))=D;
                
             end
                
                
         
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
            
             
        end
       
    end
end


        for j=y1:y2
            for i=1:nx
        
        
                if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                    continue
                
       
                elseif i==1
                       GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                       CM(f(i,j),f(i,j))=B;
                       CM(f(i,j),f(i+1,j))=2*D;
                       CM(f(i,j),f(i,j+1))=D;
                       CM(f(i,j),f(i,j-1))=D;
                   
                
                elseif i==x1-1
                       GM(f(i,j))=Tm(i,j)-D*Tin;
                       CM(f(i,j),f(i,j))=A;
                       CM(f(i,j),f(i-1,j))=D;
                       CM(f(i,j),f(i,j+1))=D;
                       CM(f(i,j),f(i,j-1))=D;
                    
                    
                
                
               
                elseif i==x2+1
                        CM(f(i,j),f(i,j))=A;
                        CM(f(i,j),f(i,j-1))=D;
                        CM(f(i,j),f(i,j+1))=D;
                        CM(f(i,j),f(i+1,j))=D;
                        GM(f(i,j))=Tm(i,j)-D*Tin;
                        
                    
                elseif i==nx
                        CM(f(i,j),f(i,j))=B;
                        CM(f(i,j),f(i,j-1))=D;
                        CM(f(i,j),f(i-1,j))=2*D;
                        GM(f(i,j))=Tm(i,j)+2*Bi*Fo*Ta;
                        CM(f(i,j),f(i,j+1))=D;  
                         
                         
                else
                        GM(f(i,j))=Tm(i,j);
                        CM(f(i,j),f(i,j))=A;
                        CM(f(i,j),f(i-1,j))=D;
                        CM(f(i,j),f(i,j+1))=D;
                        CM(f(i,j),f(i+1,j))=D;
                        CM(f(i,j),f(i,j-1))=D; 
                   
                   
                end
                        
                    
            end
                
        end
     

   


for j=y2+1:ny
    for i=1:nx
        
        if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                continue
   
        elseif j==y2+1
            
            if i==1
                    GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                    CM(f(i,j),f(i,j))=B;
                    CM(f(i,j),f(i+1,j))=2*D;
                    CM(f(i,j),f(i,j+1))=D;
                    CM(f(i,j),f(i,j-1))=D;
            elseif i==nx
                    GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                    CM(f(i,j),f(i,j))=B;
                    CM(f(i,j),f(i-1,j))=2*D;
                    CM(f(i,j),f(i,j+1))=D;
                    CM(f(i,j),f(i,j-1))=D;
            
            elseif i>=x1 && i<=x2
                GM(f(i,j))=Tm(i,j)-D*Tin;
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i,j+1))=D;
            else
                GM(f(i,j))=Tm(i,j);
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i,j-1))=D; 
            end
        elseif j==ny
    
            if i==1
                GM(f(i,j))=Tm(i,j)+4*Bi*Fo*Ta;
                CM(f(i,j),f(i,j))=C;
                CM(f(i,j),f(i,j-1))=2*D;
                CM(f(i,j),f(i+1,j))=2*D;

            elseif i==nx
                GM(f(i,j))=Tm(i,j)+4*Bi*Fo*Ta;
                CM(f(i,j),f(i,j))=C;
                CM(f(i,j),f(i,j-1))=2*D;
                CM(f(i,j),f(i-1,j))=2*D;
            else
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                CM(f(i,j),f(i,j))=B;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j-1))=2*D;
            end
        else
            if i==1
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                    CM(f(i,j),f(i,j))=B;
                    CM(f(i,j),f(i+1,j))=2*D;
                    CM(f(i,j),f(i,j+1))=D;
                    CM(f(i,j),f(i,j-1))=D;
            elseif i==nx
                GM(f(i,j))=Tm(i,j)+2*Fo*Bi*Ta;
                    CM(f(i,j),f(i,j))=B;
                    CM(f(i,j),f(i-1,j))=2*D;
                    CM(f(i,j),f(i,j+1))=D;
                    CM(f(i,j),f(i,j-1))=D;
            else
                GM(f(i,j))=Tm(i,j);
                CM(f(i,j),f(i,j))=A;
                CM(f(i,j),f(i-1,j))=D;
                CM(f(i,j),f(i,j+1))=D;
                CM(f(i,j),f(i+1,j))=D;
                CM(f(i,j),f(i,j-1))=D;
            end
        
    end
    end

end

ANS=(CM)\GM;                                                   
clear i j
    for j=1:ny
        for i=1:nx
            if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                Tm(i,j)=T(i,j);
            else
           Tm(i,j)= ANS(f(i,j));
            end
        end 
            
    end
end



