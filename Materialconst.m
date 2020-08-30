function Csmat = Materialconst(lamda,mu,kron_delta)
Csmat = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    Csmat(i,j,k,l) = Csmat(i,j,k,l) + lamda*kron_delta(i,j)*kron_delta(k,l)...
                        + mu*(kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k));
                end
            end
        end
    end

end