
%CFD Project Number 1
%Unsteady heat transfer in a square channel with convection in outside and
                                        %constant temperature inside
%the problem is solved by 3 Methods: FTCS - PSOR - LASSONEN
%different values of time grid has been set and plotted
%Dt is optional - If it is changed Pay attention to plots and stability



%%
close all
clear 
clc
%Student No. :401742274

a0=1+0.02*74;                      %Channel inner width   
ai=2*a0;                           %Channel outer width

W=abs(ai-a0)/2;                    %Channel Wall
w=3;                               %number of partitions in Channel Wall
dx=W/w;                            %delta x
dy=dx;                             %delta y

alpha=17.7e-6;                     %Diffusion Coefficient
h=50;                              %Convection Coefficient
k=60.5;                            %Conductivity

tf=3600*20;                        %Final Time (if needed-not used here)
x=0:dx:ai;                         %x vector
y=0:dy:ai;                         %y vector
nx=length(x);                      %Number of nodes in x direction
ny=length(y);                      %Number of nodes in x direction
Bi=h*dx/k;                         %finite-difference form of the Biot number 

T0=300;                            %Initial Temperature
T=ones(nx,ny)*T0;                  %Initializing Temperature
Tin=400;                           %inner Temperature
Ta=300;                            %outer Fluid Temperature

x1=find(x==W); x2=find(x==(ai-W)); %Position Of Inner Corners
y1=find(x==W);y2=find(x==(ai-W));

%% Initializing Temperature
for i=1:nx
    for j=1:ny
        if i>=x1 && i<=x2 && j==y1
            T(i,j)=Tin;
        elseif i>=x1 && i<=x2 && j==y2
            T(i,j)=Tin;
        elseif j>=y1 && j<=y2 && i==x1
            T(i,j)=Tin;
        elseif j>=y1 && j<=y2 && i==x2
            T(i,j)=Tin;
        elseif (i>x1 && i<x2) && (j>y1 && j<y2)
            T(i,j)=nan;
        end
        
    end
end
%% PSOR Method
[T_PSOR, it_PSOR] = PSOR(T,nx,ny,x1,x2,y1,y2,Ta,Bi);
%%
q=1;                                      %counter
Dt=[50 100 120 130 150 190  ];                      %delta t
T_FTCS=cell(1,length(Dt));                %%Initializing a cell to reserve Temperatures With Diiferent delta ta
T_LASSONEN=cell(1,length(Dt));

%% FTCS and LASSONEN Method solution
for dt=Dt
t=0:dt:tf;                                     %time vector
nt=length(t);                                  %number of time nodes

Fo=alpha*dt/(dx)^2;                            %finite-difference form of the Fourier number
St1=Fo*(1+Bi);                                 %Stability Condtion
St2=Fo*(2+Bi);                                 %Stability Condtion

while Fo>1/4 || St1>1/4 || St2>1/2             %Checking Stability
    disp('Not Stable')
    dt=input('Insert new dt:\n');
    Fo=alpha*dt/(dx)^2;
    St1=Fo*(1+Bi);
    St2=Fo*(2+Bi);
    t=0:dt:tf;
    nt=length(t);
end


eps=10^-4;                                     %epsilon
T_FC_n=T;                                      %new FTCS Temp. matrix      

E=1000;                                        %initial diffrence of two adjacent Temp. Matrixes
it_f=0;                                        %number of iteration Until Steady State is rechead

while E>eps
    T_FC=T_FC_n;
    T_FC_n = FTCS(T_FC,nx,ny,x1,x2,y1,y2,Fo,Bi,Ta);      %FTCS solution
    E=max(max(abs (T_FC_n-T_FC  )));                     %Maximum diffrence of two adjacent Temp. Matrixes
    it_f=it_f+1;
end
T_FC=T_FC_n;

T_FTCS{q}=T_FC;                                           %saving Temp. in a cell


FT_F=(it_f-1)*dt;                                         %time when steady state is reached in second
FT_FTCS_hour(q)=FT_F/3600;        


%{                                                 
for i=1:nx                                           plotting Mesh
    for j=1:ny
        if isnan(T(i,j))
            MESH(i,j)=nan;
        else 
            MESH(i,j)=1;
            plot(i,j,'r.')
            hold on
        end
    end
end
%}                     


Tm_N=T;
E2=100;                                                 %initial diffrence of two adjacent Temp. Matrixes
it_l=0;                                                 %number of iteration Until Steady State is rechead
while E2>eps
    Tm=Tm_N;                         
    Tm_N = Lassonnen(Tm,nx,ny,x1,x2,y1,y2,Ta,Tin,w,Fo,Bi);  %LASSONEN solution
    E2=max(max(abs(Tm_N-Tm)));                              %Maximum diffrence of two adjacent Temp. Matrixes
    it_l=it_l+1; 
end
Tm=Tm_N;
T_LASSONEN{q}=Tm;                                                %saving Temp. in a cell
FT_L=(it_l-1)*dt;                                                %time when steady state is reached in second
FT_LASSONEN_hour(q)=FT_L/3600;
Error_Lassonen(q)=max(max(abs((Tm-T_PSOR)./T_PSOR)))*100;       %maximum error for FTCS with respect to PSOR 
Error_FTCS(q)=max(max(abs((T_FC-T_PSOR)./T_PSOR)))*100;         %maximum error for LASSONEN with respect to PSOR
DIFF(q)= max(max(abs((T_FC-Tm))));

disp(['Time until steady state is reached in FTCS Method in hour: ',num2str(FT_FTCS_hour(q))]);
disp(['Time until steady state is reached in Lassonen Method in hour: ',num2str(FT_LASSONEN_hour(q))]);
q=q+1;
end

FT_FTCS_hour=FT_FTCS_hour';
FT_LASSONEN_hour=FT_LASSONEN_hour';
Error_Lassonen=Error_Lassonen';
Error_FTCS=Error_FTCS';
DIFF=DIFF';

%% PLOTTING
[X,Y]=meshgrid(x,y);

figure(1)
for i=1:length(Dt)
    subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i)
    contourf(X,Y,cell2mat(T_FTCS(i)))
    title(['Temperature Disturbution(℃), FTCS, Δt= ', num2str(Dt(i)),' s'])
    ylabel('Height (m)');xlabel('Width (m)');
    colorbar
    set(gca,'fontsize',11);
    hold on
end

figure(2)
for i=1:length(Dt)
    subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i);
    TT=cell2mat(T_FTCS(i));
    plot(x,TT(round(ny/2),:),'Linewidth',2.5)
    title(['Temperatue at y=',num2str(y(round(ny/2))),', Δt=',num2str(Dt(i)),' s']) 
    xlabel('Width (m)');ylabel('Temperature (℃)');xlim([0 ai]);
    set(gca,'fontsize',14);
    hold on
end

figure (3)
for i=1:length(Dt)
    subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i)
    contourf(X,Y,cell2mat(T_LASSONEN(i)))
    title(['Temperature Disturbution(℃), LASSONEN, Δt= ', num2str(Dt(i)),' s'])
    ylabel('Height(m)');xlabel('Width (m)')
    set(gca,'fontsize',11);
    colorbar
    hold on
end

figure(4)
for i=1:length(Dt)
    subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i);
    TT2=cell2mat(T_LASSONEN(i));
    plot(x,TT2(round(ny/2),:),'Linewidth',2.5)
    title(['Temperatue at y=',num2str(y(round(ny/2))),', Δt=',num2str(Dt(i)),' s']) 
    xlabel('Width (m)');ylabel('Temperature (℃)'); xlim([0 ai]);
    set(gca,'fontsize',14);
    hold on
end

figure(5)
 contourf(X,Y,T_PSOR)
    title(['Temperature Disturbution(℃), PSOR'])
    ylabel('Height(m)');xlabel('Width (m)')
    set(gca,'fontsize',18);colorbar;

    
    %% One Fourth of the channel            
    
   
    

    i=1;
    T4=cell(1,length(Dt));
    T_one_fourth=cell(1,length(Dt));
    for dt=Dt
        [T4{i}, T_one_fourth{i}, FT_onefourth(i), X2, Y2,m]=one_fourth(a0,ai,w,alpha,h,k,T0,Tin,Ta,eps,dt);
        Error_onefourth(i)=max(max(abs((T4{i} - cell2mat(T_FTCS(i)))./cell2mat(T_FTCS(i)))))*100;
        i=i+1;
        
    end

    
  for i=1:length(Dt)
    
        figure(6)
        subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i)
        contourf(X,Y,T4{i})
        title(['Temperature Disturbution(℃), FTCS from One Fourth of channel, Δt= ', num2str(Dt(i)),' s'])
        ylabel('Height (m)');xlabel('Width (m)');
        colorbar
        set(gca,'fontsize',11);
        hold on
        figure(7)
        subplot(ceil(length(Dt)/2),ceil(length(Dt)/2),i)
        contourf(X2,Y2,T_one_fourth{i})
        title(['1/4 Channel Temperature Disturbution(℃), FTCS, Δt= ', num2str(Dt(i)),' s'])
        xlabel('Width (m)');ylabel('Heigth (m)');
        hold on
  end
  
  