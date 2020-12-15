%  ======================================================================
%              Two-dimensional Ising model 
%    Phase Transition from paramagnetism to ferromagnetism 
%  ========================================================================

% In this model we take value Kb=1 and J=1
clear all;
close all;

n = 25;  %No of lattice points.
m= ones(n); % forming a intital lattice of all spin up.

iter = 200; % no. of iterations.

% Boundary condition for positive and negative x,y axis
px(1:n)=(1:n)+1; px(n)=1;
nx(1:n)=(1:n)-1; nx(1)=n;
py(1:n)=(1:n)+1; py(n)=1;
ny(1:n)=(1:n)-1; ny(1)=n;

%forming empty arrays for different values. 
eng=[];
neigh=[];
avg_eng=[];
Cv=[];
sus=[];

%set initial and final temperature
in_temp = 1;
en_temp = 5;

%initiating the loop with temperature difference 0.01 
for t=en_temp:-0.01:in_temp
    for hh =1:iter           %loop for no. of iteration perform
        for i=1:n            %loop for lattice site i
            for j=1:n        %loop for lattice site j
                neigh(i,j) = m(px(i),j)+m(nx(i),j)+m(i,py(j))+m(i,ny(j));  %calculating the sum of all the neighbour of lattice site (i,j).
                dE = 2*neigh(i,j)*m(i,j);                                  %calculating the change in energy.
                
%                 =======================
                   %metropolis algorithm
%                 =======================
                
                r=rand();                           % generating a random number.
                if dE <0 || exp(-dE/t)>r            % checking if change in energy is less than zero or random no. generated is less than boltzmann factor.
                m(i,j)=-m(i,j);                     % if above statement is true then flip the spins.
              
                end
            end
        end
    end
    figure(1);
    imagesc(m)
    pause(0.00000000001)
eng=-0.5*(m.*neigh);                                %cal the energy
avg_eng = mean(eng,'all');                          %cal the average energy
avg_mag = sum(m,'all')/(n^2);                       %cal the average magnetisation
mag = sum(m.*m,'all')/(n^2);                        %cal the magnetisation 
Cv = (mean((eng.*eng),'all')-(avg_eng.*avg_eng))/(t*t);      %cal the specific heat
sus=(mean((mag.*mag),'all')-(avg_mag.*avg_mag))/(t);         %cal the magnetic susceptibility

%Plotting the Average Energy,Average Magnetisation,Specific Heat,Magnetic
%Susceptibilty w.r.t Temperature
figure(2);
subplot(2,2,1); 
plot(t,avg_eng,'r*');title('Average energy Vs Temperature');xlabel('Temperature');ylabel('Average energy');
hold on;
subplot(2,2,2);plot(t,avg_mag,'g*');title('Average Magnetisation Vs Temperature');xlabel('Temperature');ylabel('Average Magnetisation');
hold on;
subplot(2,2,3);plot(t,Cv,'b*');title('Specific Heat Vs Temperature');xlabel('Temperature');ylabel('Specific Heat');
hold on;
subplot(2,2,4);plot(t,sus,'k*');title('Susceptibility Vs Temperature');xlabel('Temperature');ylabel('Magnetic Susceptibility');
hold on;


end

