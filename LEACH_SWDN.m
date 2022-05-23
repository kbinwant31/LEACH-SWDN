%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%                     LEACH-SWDN Implementation                        %
%                                                                      %                                  
%     "A clustering algorithm based on energy information              %
%                      and cluster heads                               %
%            expectation for wireless sensor networks "                %                                                             
%                                                                      %  
%                                                                      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBMITTED BY-                                                        %
%                02001012019: MITALI LAROIA                            %
%                03001012019: BINWANT KAUR                             %
%                03801012019: ANJALI SINGH                             %
%                04001012019: AASTHA CHAUDHARY                         %
%                     (B.Tech CSE-1)                                   %
%                                                                      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=1.75*ym;

%Number of Nodes in the field
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0.05;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.05;
%\alpha
a=0.1;

%maximum number of rounds
%rmax=1500
rmax=1500
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+');
        hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
 
        
%First Iteration
figure(1);

%counter for CHs per round
rcountCHs=0;
%counter for CHs
countCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

enrgy_res=zeros(n,rmax)+Eo;
energy_moy=zeros(n,rmax)+Eo;
Energy_disp=zeros(n,rmax);

%last round
last=rmax;

for r=1:1:rmax
    r

  %Operation for epoch
  if(mod(r-1, round(1/p) )==0)
    for i=1:1:n
      S(i).G=0;
      S(i).cl=0;
    end
  end

  hold off;


 %Number of dead nodes
 dead=0;
 %Number of dead Advanced Nodes
 dead_a=0;
 %Number of dead Normal Nodes
 dead_n=0;
 %Alive node
 alive=0;

 %counter for bit transmitted to Bases Station and to Cluster Heads
 packets_TO_BS=0;
 packets_TO_CH=0;
 %counter for bit transmitted to Bases Station and to Cluster Heads 
 %per round
 PACKETS_TO_CH(r)=0;
 PACKETS_TO_BS(r)=0;



 figure(1);

 for i=1:1:n
    
   %checking if there is a dead node
   if (S(i).E<=0)
     plot(S(i).xd,S(i).yd,'black .');
     dead=dead+1;
     if(S(i).ENERGY==1)
       dead_a=dead_a+1;
     end
     if(S(i).ENERGY==0)
       dead_n=dead_n+1;
     end
     hold on;    
   end
   if S(i).E>0
     S(i).type='N';
     if (S(i).ENERGY==0)  
       plot(S(i).xd,S(i).yd,'o');
     end
     if (S(i).ENERGY==1)  
       plot(S(i).xd,S(i).yd,'+');
     end
     hold on;
   end
 end
 plot(S(n+1).xd,S(n+1).yd,'red x');


 STATISTICS(r+1).DEAD=dead;
 DEAD(r)=dead;
 DEAD_N(r)=dead_n;
 DEAD_A(r)=dead_a;
 alive=n-dead;
 ALIVE_NODE(r)=alive;

 if alive<10
     last=r;
     break;
 end

  %When the first node dies
  if (dead==1)
    if(flag_first_dead==0)
      first_dead=r;
      flag_first_dead=1;
    end
  end

 countCHs=0;
 cluster=1;


 for i=1:1:n
   if(S(i).E>0)

     %average energy
     energy_moy(i,r)=(1/r)*sum(enrgy_res(i,1:(r)));
     % the rundom number is generated at the intervall [0 ,(Eaverage_energy(i)/Eo) ]
     % generate a rundom number
     temp_rand= ((energy_moy(i,r)/Eo)* rand(1,1)); 

     if ( (S(i).G)<=0)
     
       %k = (a live node in the network * %c(number of cluster heads defined at the intiale time))
       k=ALIVE_NODE(r)*p;
       %Election of Cluster Heads
       if(temp_rand<= ( k*(S(i).E/Eo) / (n-k*mod(r,round(n/k))) ) )
         
         countCHs=countCHs+1;
         packets_TO_BS=packets_TO_BS+1;
         PACKETS_TO_BS(r)=packets_TO_BS;
            
         S(i).type='C';
         S(i).G=1/k;
         C(cluster).xd=S(i).xd;
         C(cluster).yd=S(i).yd;
         plot(S(i).xd,S(i).yd,'k*');
            
         distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
         C(cluster).distance=distance;
         C(cluster).id=i;
         X(cluster)=S(i).xd;
         Y(cluster)=S(i).yd;
         cluster=cluster+1;
            
         %Calculation of Energy dissipated
         % Size of data package is taken 4000 units in bits
         distance;
         if (distance>do)
           S(i).E=S(i).E-((ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance)); 
           enrgy_res(i,r)= S(i).E;
           Energy_disp(r)=Energy_disp(r)+((ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance));
         end
         if (distance<=do)
           S(i).E=S(i).E-((ETX+EDA)*(4000) + Efs*4000*(distance*distance)); 
           enrgy_res(i,r)= S(i).E;
           Energy_disp(r)=Energy_disp(r)+ ((ETX+EDA)*(4000) + Efs*4000*(distance*distance)); 
         end
       end  
     end
   end 
 end

 STATISTICS(r).CLUSTERHEADS=cluster-1;
 CLUSTERHS(r)=countCHs;

 %Election of Associated Cluster Head for Normal Nodes

 for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
         temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
         if ( temp<min_dis )
           min_dis=temp;
           min_dis_cluster=c;
         end
       end
       
       %Energy dissipated by associated Cluster Head
       
       if (min_dis>do)
         S(i).E=S(i).E-(ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
         enrgy_res(i,r)= S(i).E;
         Energy_disp(r)=Energy_disp(r)+ (ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
       end
       if (min_dis<=do)
         S(i).E=S(i).E-(ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
         enrgy_res(i,r)= S(i).E;
         Energy_disp(r)=Energy_disp(r)+ (ETX*(4000) + Efs*4000*( min_dis * min_dis));
       end
       %Energy dissipated
       if(min_dis>0)
         S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E-((ERX + EDA)*4000); 
         PACKETS_TO_CH(r)=n-dead-cluster+1; 
         Energy_disp(r)=Energy_disp(r)+(ETX*(4000) + Efs*4000*( min_dis * min_dis));
       end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
     end
   end
 end
 hold on;

 countCHs;
 rcountCHs=rcountCHs+countCHs;

sum_=0;
nrg=0;
for i=1:1:n
if(S(i).G<=0)
    sum_=sum_+S(i).E;
end
end

avg=sum_/n;
avg_r(r)=avg;
energy_moy(r+1)=avg;

energy_(r)=sum(enrgy_res(1:(n),r))/n;

warning('OFF');
[vx,vy]=voronoi(X(:),Y(:));
plot(X,Y,'g+',vx,vy,'m-');
title 'Wireless Sensor Network';
hold on;
voronoi(X,Y);
axis([10 xm 0 ym]);

end

%%%%%%%%%%%%%%%%%%%% Plotting Simulation Results %%%%%%%%%%%%%%%%%%%%

    % "Operating Nodes per Round" %
    figure(2)
    plot(1:r,ALIVE_NODE(1:r),'-r','Linewidth',2);
    title ({'LEACH-SWDN'; 'Network lifetime';})
    xlabel 'Time';
    ylabel 'No of alive nodes';
    hold on;
    

    figure(3)
    plot(1:last,cumsum(Energy_disp(1:last)),'-r','Linewidth',2);
    title ({'LEACH-SWDN'; 'Total Energy Consumption of Nodes';})
    xlabel 'Time';
    ylabel 'Energy (J)';
    hold on;

    figure(4)
    plot(1:last,CLUSTERHS(1:last),'-*','MarkerIndices',1:last);
    title ({'LEACH-SWDN'; 'Cluster Heads per Round';})
    xlabel 'Round';
    ylabel 'No of Cluster Heads';
    hold on;

    figure(5)
    plot(1:last,cumsum(PACKETS_TO_BS(1:last)),'-r','Linewidth',2);
    title ({'LEACH-SWDN'; 'Total Packets received at Base Station';})
    xlabel 'Time';
    ylabel 'No of Packets received';
    hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 							  %
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round					  %
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round                     %
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round                      %
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round      %
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round   %
%  first_dead: the round where the first node died                                    %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%