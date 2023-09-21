function A = dmd_recon_r0(Phi, delta, omega, dt, K)
    time_dynamics = zeros( size(Phi,2), K );
    t = dt*linspace(0, K-1, K);
    for k=K:-1:1
       time_dynamics(:,k) = exp( (delta+1i*omega)*(t(k)-t(1)) );
    end
    A = real( Phi * time_dynamics ) ; %real porque devuleve x + 0*i

end