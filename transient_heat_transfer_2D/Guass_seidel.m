function [T_PSOR,it_PSOR] = Guass_seidel(T,nx,ny,x1,x2,y1,y2,Ta,Bi)


eps=0.001; 
error2=10;
T_PSOR=T;             %gauss-seidel temperature matrix
it_PSOR=0;            %number of iterations until convergence
w=1.5;
while error2>eps

    T_old_PSOR=T_PSOR;       %T_GS at k-1 iteration
    
    for i=1:nx
        for j=1:ny
           
            
            if (i>=x1 && i<=x2) && (j>=y1 && j<=y2)
                continue
          
            elseif i==1 && j==1
                T_PSOR(i,j)= w/(2*Bi+2) * (T_old_PSOR(i+1,j)+T_old_PSOR(i,j+1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            elseif i==1 && j>1 && j<ny
                T_PSOR(i,j)=w/(2*Bi+4)*(2*T_old_PSOR(i+1,j)+T_old_PSOR(i,j+1)+T_PSOR(i,j-1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
                
            elseif j==ny && i==1
                T_PSOR(i,j)=w/(2*Bi+2) *(T_old_PSOR(i+1,j)+T_PSOR(i,j-1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            elseif i==nx && j==1
                T_PSOR(i,j)=w/(2*Bi+2)*(T_PSOR(i-1,j)+T_old_PSOR(i,j+1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            elseif i==nx && j==ny
                T_PSOR(i,j)=w/(2*Bi+2)*(T_PSOR(i-1,j)+T_PSOR(i,j-1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
                
                
            elseif j==1 && i>1 && i<nx
                T_PSOR(i,j)=w/(2*Bi+4)*(2*T_old_PSOR(i,j+1)+T_old_PSOR(i+1,j)+T_PSOR(i-1,j)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            elseif j==ny && i>1 && i<nx
                T_PSOR(i,j)=w/(2*Bi+4)*(2*T_PSOR(i,j-1)+T_old_PSOR(i+1,j)+T_PSOR(i-1,j)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            
            elseif i==nx && j>1 && j<ny
                T_PSOR(i,j)=w/(2*Bi+4)*(2*T_PSOR(i-1,j)+T_old_PSOR(i,j+1)+T_PSOR(i,j-1)+2*Bi*Ta)+(1-w)*T_old_PSOR(i,j);
            else
                T_PSOR(i,j)=w/4 *(T_old_PSOR(i+1,j)+T_PSOR(i-1,j)+T_PSOR(i,j-1)+T_old_PSOR(i,j+1))+(1-w)*T_old_PSOR(i,j);
            end
        end
    end
    
    error2=max(max(abs(T_PSOR-T_old_PSOR)));
    it_PSOR=it_PSOR+1;
end


end

