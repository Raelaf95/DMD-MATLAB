function [lambda, Phi, Amplitude, delta, omega, f] = dmd_d_r0(X, d, eps, epsilon_spectra, dt)
%     [Lambda, Phi, Amplitude, delta, omega, f] = dmd_d_r0(X, 3, 1E-3, 1E-3, 1/24)
    %% STEP 1: tuncate original data
    [U_Xr, S_Xr, V_Xr] = trunc_svd(X, eps);
%     Xr = U_Xr * S_Xr * V_Xr';
    [J,K]=size(X);
%     disp( ["err. trunc. matriz original " norm(X-Xr)] )
    %% STEP 2: Create reduced snapshots matrix
    Xhat = S_Xr * V_Xr';
    %% STEP 3 Create the modified snapshot matrix
    [N,~]=size(Xhat);
    Xtilde=zeros(d*N,K-d+1);
    for iii=1:1:d
        Xtilde((iii-1)*N+1:iii*N,:) = Xhat(:,iii:iii+K-d);
    end
    %% Step 4 Dimension reduction of the modified snapshot matrix
    [U_Xtilder, S_Xtilder, V_Xtilder] = trunc_svd(Xtilde, eps);
    Xtildehat = S_Xtilder * V_Xtilder';
    %% Step 5 Reduced modified snapshot matrix
    X1 = Xtildehat(:,1:end-1);
    X2 = Xtildehat(:,2:end);
    %% Step 6 Hypothesis: X2 = R * X1 working with reduced spaces, I can 
    %% Step 7 calculate Rtilde = U' * R * U = U' * X2 * V *S^-1
    [Ur, Sr, Vr] = svd(X1,'econ');
%     R = Ur' * X2 * Vr / Sr;
%     [Q, Lambda] = eig(Rtilde);    
    % %%% Calculo R
    Rtilde = X2 * Vr / Sr * Ur';
    [tildeQ, Lambda] = eig(Rtilde);   
    %% Step 8 calculate omega and delta
    lambda = diag(Lambda);
    delta = real(log(lambda))/dt;
    omega  = imag(log(lambda))/dt;
    %% Step 9 calculate amplitudes
    Q = U_Xtilder * tildeQ;
    Q = Q( (d-1)*N+1:d*N,: ); %autovectores del último ciclo
    [NN,MM]=size(Q);
    for j = 1:1:MM %ortonormal
        Q(:,j)= Q(:,j)/norm(Q(:,j));
    end
    mat_L = zeros(NN*K,MM);
    mat_L(1:NN,:) = Q*eye(MM);
    for kk = 2:1:K
        mat_L((kk-1)*NN+1:kk*NN,:) =  mat_L((kk-2)*NN+1:(kk-1)*NN,:)*Lambda;
    end
    vec_b = Xhat(:);
    %%% L * a = b -> a = psinv(L) * b
    [u,s,v] = svd(mat_L,'econ');
    psimat_L = v .* (1./diag(s))' * u' ;
    vec_a = psimat_L * vec_b;
    mat_u = reshape(vec_a,[1, numel(vec_a)]).*Q;%
    
    %%% convierto vec_a en ai*ui
    Phi = U_Xr * mat_u;%reshape(vec_a,[1, numel(vec_a)]).*mat_u;
    % hatPhi = Ur * ( reshape(vec_a ,[1, numel( vec_a )]) .* Q * Ur' );
    for j= MM:-1:1
        Amplitude(j)=norm( Phi(:,j) ) / sqrt(J);
    end
    %% Step 10 Ordenar y truncar modos
    [Amplitude,Idx_A] = sort(Amplitude,'descend');
    Phi = Phi(:,Idx_A);
    delta = delta(Idx_A);
    omega  = omega(Idx_A);
    lambda = lambda(Idx_A);
%     mat_u = mat_u(:,Idx_A);
    %%% Truncar modos
    NmodosDMD = find(Amplitude/Amplitude(1) < (epsilon_spectra), 1, 'first');
    if isempty(NmodosDMD), NmodosDMD = MM;
    else
        NmodosDMD = NmodosDMD - 1;
    end
    disp([ 'Number of DMD modes: ' num2str(NmodosDMD) ])
    Phi = Phi(:,1:NmodosDMD);
    delta=delta(1:NmodosDMD);
    omega=omega(1:NmodosDMD);
    Amplitude=Amplitude(1:NmodosDMD)';
%     mat_u = mat_u(:,1:NmodosDMD)';
    f = omega/2/pi;
    disp(array2table([(1:NmodosDMD)',delta,omega,Amplitude],...
        'VariableNames',{'Mode', 'delta', 'omega', 'Amplitude'}))
    %% Step 11 DMD modes
%     Phi = zeros(J,NmodosDMD);
%     for m=1:1:NmodosDMD
% %         NormMode = norm(U_Xr*hatPhi(:,m),2)/sqrt(J);
%         Phi(:,m) = hatPhi(:,m)/NormMode;
%         Phi(:,m) = Amplitude(j)*hatPhi(:,m);
%         Phi(:,m) = J * hatPhi(:,m) / norm(hatPhi(:,m));
%     end
end