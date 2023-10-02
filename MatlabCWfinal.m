%MATLAB coursework - Aeroelasticity 
%Code which runs the simulation of an airfoil for two different speeds
%And plots pitch and heave for different times 

%Clearing
clear
clc 

%Defining variables
c = 1.5;              % m 
m = 10;               % kg/m
dt = 0.01;            % s
dhdt = 0;             % m/s
dadt = 0;             % rad/s
w1 = 0;                
w2 = 0;               
w3 = 0;               
w4 = 0;               
h_0 = 0;              % m
xc = 0.2;             
xf = 0.22;            
I_a = 100;            % kgm
K_h = 400;            % N/m
rho = 1.225;          % kg/m^3 
K_alpha = 500;        % Nm/rad
alpha_0 = deg2rad(5); % rad

x = [0; 0; 0; alpha_0; 0; 0; 0; 0];

%Defining Velocities, v
%Velocity 1 = 9.5 m/s, Velocity 2 = 12 m/s 
v = 12;

%Defining Matrices M and K
M = getM(m,rho,c,xc,xf,I_a);
K = getK(v,K_h,K_alpha,c,xf,rho);
%Inversing Matrix M and multiplying by K
invMK = M\K;

%Defining time frame 
time = 0: 0.01: 45;
l = length(time); 

%Loading Data and defining positions of airfoil
load('NACA0012 (1).mat','xairf','yairf')
xpos = xairf';
ypos = yairf';

%Creating a figure
figure
line = animatedline('color','m','Linewidth',2);
varv = sprintf('%.2f',v);
%Adding title and setting axis 
title(['Aeroelastic simulation with v = ' varv, 'm/s'], 'FontSize', 25)
axis([-0.4 1.4 -0.5 0.5])
%Labelling the axes of plots and plotting airfoil
xlabel('Length (m)', 'FontSize', 16)
ylabel('Variation', 'FontSize', 16)
grid on
set(gca, 'Fontsize', 16); 
hold on 

%Creating video (vid) and saving it in the correct format 
vid = VideoWriter(['Airfoil' varv, '.avi'], 'Uncompressed AVI');
vid.FrameRate = 100;
open(vid)

%Creating pitch and heave arrays 
pitch = (zeros(l,1))';
heave = (zeros(l,1))'; 

%Defining index to be used in loop 
ii = 1;

%Use of Runge Kutta condition to find k_1 until k_4 and new x 
for i = time 
    k_1 = invMK*x*dt;
    k_2 = invMK*(x+0.5*k_1)*dt;
    k_3 = invMK*(x+0.5*k_2)*dt;
    k_4 = invMK*(x+k_3)*dt;
    x = x + (k_1 + 2*k_2 + 2*k_3 + k_4)/6;

    pitch(ii) = x(4); % Storing the new values for pitch and heave
    heave(ii) = x(3);
    ii = ii + 1;

    %New positions
    xn = (xpos - xf)*c*cos(x(4)) + ypos*c*sin(x(4));
    yn = -x(3) - (xpos - xf)*c*sin(x(4)) + ypos*c*cos(x(4));

    %Creating the video
    clearpoints(line)
    addpoints(line,xn,yn)
    drawnow
    pause(0.001)
    legend(sprintf('%.2fs',i))
    frame = getframe(gcf);
    writeVideo(vid,frame)

end 

%Closing the video
close(vid)
    
%Pitch-time plot 
figure
subplot(2,1,1);
plot(time,pitch,'Color', 'k', 'LineWidth',2);
axis([0 45 -0.15 0.15])
%Labelling
xlabel('Time (s)', 'FontSize', 20)
ylabel('Pitch (°)', 'FontSize', 20)
set(gca, 'Fontsize', 16); 
title(['How pitch varies with time with v = ' num2str(v), 'm/s'], 'FontSize',18)

%Heave-time plot 
subplot(2,1,2);
plot(time,heave,'Color','k','LineWidth',2);
axis([0 45 -0.15 0.15])
%Labelling
xlabel('Time (s)', 'FontSize', 20)
ylabel('Heave (°)', 'FontSize', 20)
set(gca, 'Fontsize', 16); 
title(['How heave varies with time with v = ' num2str(v), 'm/s'], 'FontSize',18)

%Changing axis conditions depending on value for velocity 
if v == str2double('12.0')
    axis([0 45 -0.4 0.4])
end 

%Saving plot 
saveas(gcf,['Pitch_and_Heave_plot_for ' varv '.png'],'png');

%End of script