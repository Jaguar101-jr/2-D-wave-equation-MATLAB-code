%% This code is written by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED, B.Sc. UofK-Sudan,
%% now at UPC (China), August, 2021 as part of his master thesis's project.
%% A Finite Difference Method program to solve/propagate the 2D constant-density,
%% wave equation on a Cartesian mesh, wth 2D, 3D plots and a movie.
%% The 2D Acoustic Wave Equation: Wtt = c^2 (Wxx + Wyy) + f(t).
%% Please cite this code as: Hussein Abduelhaleim Hussein Muhammed (2022), Least-
%% Squares Reverse-Time Migration in Pseudodepth Domain and Its Applications, 
%% China University of Petroleum (East China), Master Thesis,
%% School of Geosciences, Dept. of Geophysics Library press.
%% Thanks is to Dr. Haroon Stephen (UNLV-Las Vegas, USA) and Dr. Gao Liping (UPC
%% Qingdao-China).

clear;

%% Prepare the movie file.
    vidObj = VideoWriter('twoDwave.avi');
    open(vidObj);
 
%% Creating the Domain
% Creating Space
Lx=10;
Ly=10;
dx=0.1;
dy=dx;
nx=fix(Lx/dx);
ny=fix(Ly/dy);
x=linspace(0, Lx, nx);
y=linspace(0, Ly, ny);

%Wavefield total propagation time
T=20;

%% Field variables
% Variables
wn=zeros(nx,ny);
wnm1=wn; % w at time n-1
wnp1=wn; % w at time n+1

% Parameters for stability conditions
% 0.5>=CFL=<0.6
CFL=0.55; % CFL = c.dt/dx
c=1;
dt=CFL*dx/c;

%% No Initial Conditions

%Time stepping Loop
t=0;

while(t < T)
    
%% Activate any boundary conditions
   % Reflecting Boundary Conditions
    %wn(:,[1 end])=0;
    %wn([1 end],:)=0;
   % Mur's Absorbing Boundary Conditions
    wnp1(1,:)=wn(2,:) +((CFL-1)/(CFL+1))*(wnp1(2,:)-wn(1,:));
    wnp1(end,:)=wn(end-1,:) + ((CFL-1)/(CFL+1)*(wnp1(end-1,:)-wn(end,:)));
    wnp1(:,1)=wn(:,2) + ((CFL-1)/(CFL+1)*(wnp1(:,2)-wn(:,1)));
    wnp1(:,end)=wn(:,end-1) + ((CFL-1)/(CFL+1)*(wnp1(:,end-1)-wn(:,end)));
   
    
    %Time Solution
    t=t+dt;
    wnm1=wn; wn=wnp1; % Save current and previous arrays
    
    % Source Wavelet
    wn(50,50)=dt^2*20*sin(30*pi*t/20);
    %Finite Difference Solution loop
    for i=2:nx-1,  for j=2:ny-1
        wnp1(i,j) = 2*wn(i,j)- wnm1(i,j) ...
                    + CFL^2 * (wn(i+1,j) + wn(i,j+1) - 4*wn(i,j) + wn(i-1,j) + wn(i,j-1));
    end, end

    % Check Convergence
    
    % Visualize at selected steps
    % Creating a wavefield movie
    clf;
    %subplot(2,1,1);
    %imagesc(x, y, wn'); colorbar; caxis([-0.02 0.02])
    %title(sprintf('t = %.2f' , t));
    %subplot(2,1,2);
    mesh(x, y, wn'); colorbar; caxis([-0.02 0.02])
    title(sprintf('XYZ: fourth_axis_time_elapsed = %.2f' , t));
    axis([0 Lx 0 Ly -0.02 0.02]);
    shg; pause(0.01);
    
       % Write each frame to the file.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);
