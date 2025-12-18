% Number of gridpoints and domain size
N = 200; L = 10;
% position of kink
a=0;


% Spatial grid and grid spacing
x = linspace(-L,L,N);
dx = x(2)-x(1);

v=0;
gamma = 1/sqrt((1-v^2));

% Initial condition
p0 = tanh(gamma*(2*x-a)) ;


% Initialize the spatial parameters for simulation
T = 51; dt = 0.01; time = 0:dt:T;

% Preallocate the solution array
phi = zeros(N, length(time));
p1 = tanh(gamma*(2*x-a-v*dt));


phi(:, 1) = p0;
phi(:, 2) = p1; % Set initial condition


% Time-stepping loop for the simulation
for n = 2:length(time)-1
    % Update the solution using a numerical scheme (e.g., explicit Euler method)
   k=2:length(phi(:,1))-1;
   phi(1,n+1) = phi(1,n) + (dt/dx)*(phi(2,n) - phi(1,n));
   phi(end,n+1) = phi(end,n) - (dt/dx)*(phi(end,n) - phi(end-1,n));
   k = 2:N-1;
   phi(k,n+1) = 2*phi(k,n)-phi(k,n-1)+dt^2*((phi(k+1,n)-2*phi(k,n)+phi(k-1,n))/dx^2-2*phi(k,n).*(phi(k,n).^2 - 1));
end


close all;


figure;
vid = VideoWriter('settling_kink.mp4', 'MPEG-4');
vid.FrameRate = 30;
open(vid)
saveToGif(x,phi,vid,[-L+0.1,L-0.1,-2,2],150,v,gamma,T)


figure;
plot(x,phi(:,1),LineWidth=2)
hold on;
plot(x,phi(:,round(end/2)),LineWidth=2)
plot(x,phi(:,end),LineWidth=2)
axis([-9.9,9.9,-2,2])
plot(x,tanh(gamma*x),'--'
legend('$t=0$','$t=T/2$','$t=T$','$\tanh(x)$','Interpreter','latex','Location','best','FontSize',20)
xlabel('$x$','Interpreter','latex','FontSize',24)
ylabel('$\phi$','Interpreter','latex','FontSize',24)

function saveToGif(x,phi,vid,axs,frames,v,gamma,T)
    plt1 = plot(x,phi(:,1),'LineWidth',2);
    hold on;
    plt2 = plot(x,tanh(gamma*x),'r--');
    axis(axs)
    xlabel("$x$",Interpreter='latex',FontSize=14)
    ylabel("$\phi$",Interpreter='latex',FontSize=14)
    title('Settling of Energised Kink','FontSize',20)
    for i = 1:frames
        frame= getframe(gcf);
        writeVideo(vid, frame.cdata);
        axis(axs)
        plt1.YData = phi(:,round(i*end/frames)); % Update plot data for the current frame
        plt2.YData = tanh(gamma*(x-v*round(i)*T/frames*ones(size(x)))); % Update plot data for the current frame
        drawnow; % Update the figure window
        %pause(0.1); % Pause to control the frame rate
    end
end

close(vid);
