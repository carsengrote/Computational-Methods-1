% Kuramoto–Sivashinsky equation u_t +αu_xx +βu_xxxx +u(u_x) = 0 on periodic domain [0, 2pi) solved
% with spectral method and FFT with IC u(x,0) = exp(-10 - 11cos(x))


L = 2*pi; % Length of the interval, [0,L]
N = 256; % Number of points spatially
h = L/N; % Space discretization to evaluate function at
% Fundamental frequency is 2pi/L
x = 0:h:L - h; % Discretization of interval


alpha = .1;
beta = 10^(-3);

createVideo = 1; 
if (createVideo == 1)
    writerObj = VideoWriter('testFile.mp4','MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = 20;
    open(writerObj);
end

T = 20; % Final time
k = 10^(-5); % Time step
% Initial conditions at t = 0
u = exp(-10-11*cos(x)); % 

a = fft(u); % Get the fourier coeffficients, if we have u at N points
            % then we get the a_k for -N/2 <= k <= N/2 - 1
      

uPrime = zeros(1,length(a)); % Placeholder for the derivative values
u_x = zeros(1,length(a)); % For convolution spatial derivative

%count = 0;
for n = 0: T/k - 1

    % RK4!
    k1 = Derive2(a, uPrime, u_x, alpha, beta);
    k2 = Derive2(a+ (2*k/5)*k1, uPrime, u_x, alpha, beta);
    k3 = Derive2(a+ (-3*k/20)*k1 + (3*k/4)*k2 , uPrime, u_x, alpha, beta);
    k4 = Derive2(a+ (19*k/44)*k1 -(15*k/44)*k2 + (40*k/44)*k3, uPrime, u_x, alpha, beta);
    a = a + (k/72)*(11*k1 + 25*k2 + 25*k3 + 11*k4); % Time step

    if mod(n*k,.05) == 0 % Can change the modulus number to make the video faster or slower
        
       %if count < 11
       %     plot(x,real(ifft(a)),'LineWidth',2);
       %     hold("on");
       %     count = count + 1;
       %else
       %     stop
       % end

        if createVideo == 1
            plot(x, real(ifft(a)),'LineWidth',2)
            xlabel('$$x$$','interpreter','latex');
            ylabel('$$u(x,t)$$','interpreter','latex');
            xlim([0 (L-h)]);
            ylim([-4 4]);
            set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right');
            title(['t = ', num2str(round((n*k)*100)/100) 's']);
            writeVideo(writerObj, getframe(gcf));

        end
    end

end

if createVideo == 1
     close(writerObj);
end
