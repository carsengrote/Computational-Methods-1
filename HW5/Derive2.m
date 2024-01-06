function uPrime = Derive2(u,uPrime,u_x,alpha,beta)
% v is the viscous term, v = 0 is inviscid Burgers

    N = length(u);

    for k = -floor(N/2): floor(N/2) - 1
        
        if k < 0
            neg = 1;
        else
            neg = 0;
        end

        % Diffusion term with viscosity: V*u_xx
        uPrime(N*neg+k+1) = -alpha*((1i*k)^2)*u(N*neg+k+1) - beta*((1i*k)^4)*u(N*neg+k+1);
        % Computing u_x
        u_x(N*neg+k+1) = (1i*k)*u(N*neg+k+1);
    end

    % Lets be more efficient here
    % O(NlogN) if we do it in real space!
    uReal = ifft(u);
    u_xReal = ifft(u_x);
    convR = uReal.*u_xReal; % Computing u*u_x in real space
    convF = fft(convR); % Back to Fourier space 
    
    for k = -floor(N/2): floor(N/2) - 1
        
        if k < 0
            neg = 1;
        else
            neg = 0;
        end

        % Adding u*u_x term from convolution done in real space
        uPrime(N*neg+k+1) = uPrime(N*neg+k+1) - convF(N*neg+k+1);
    end

end
