function [U, S, V] = trunc_svd(X, eps)

    [U, S, V] = svd(X,'econ');
%     ET = sum(diag(S));
    aux = find( flip(sqrt(cumsum(flip(diag(S.^2))))) / sqrt(sum(diag(S).^2) ) <= eps, 1, 'first' );
    if isempty(aux)
        aux = length(diag(S));
    end
    U = U(:,1:aux);
    S = S(1:aux,1:aux);
    V = V(:,1:aux);
    
end