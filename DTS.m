% Function to calculate Material & Geometric Nonlinearity tensor
% Calculation of 2nd PK and Sigma
function [D,T,Sigma] = DTS(F,E, lamda, mu)
Cs = zeros(3,3,3,3);
Iden = eye(3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Cs(i,j,k,l) = Cs(i,j,k,l) + lamda*Iden(i,j)*Iden(k,l) + mu*(Iden(i,k)*Iden(j,l) + Iden(i,l)*Iden(j,k));
                end
            end
        end
    end
Spk = zeros(3,3);
Din = zeros(3,3,3,3);
    for o = 1:3
        for p = 1:3
            for q = 1:3
                for r = 1:3
                    Spk(o,p) = Spk(o,p) + Cs(o,p,q,r)*E(q,r);
                    for u = 1:3
                        for v = 1:3
                            Din(o,p,q,r) = Din(o,p,q,r) + F(o,u)*F(q,v)*Cs(u,p,v,r);
                        end
                    end
                end
            end
        end
    end
D = [Din(1,1,1,1) Din(1,1,2,2) Din(1,1,1,2) Din(1,1,2,1); Din(2,2,1,1) Din(2,2,2,2) Din(2,2,1,2) Din(2,2,2,1);Din(1,2,1,1) Din(1,2,2,2) Din(1,2,1,2) Din(1,2,2,1); Din(2,1,1,1) Din(2,1,2,2) Din(2,1,1,2) Din(2,1,2,1)];

% Tin = zeros(3,3,3,3);
% Iden = eye(3,3);
% for ii = 1:3
%     for kk = 1:3
%         for jj = 1:3
%             for ll = 1:3
%                 Tin(ii,jj,kk,ll) = Tin(ii,jj,kk,ll) + Iden(ii,kk)*Spk(jj,ll);
%             end
%         end
%     end
% end
%  T = [Tin(1,1,1,1) Tin(1,1,2,2) Tin(1,1,1,2) Tin(1,1,2,1); Tin(2,2,1,1) Tin(2,2,2,2) Tin(2,2,1,2) Tin(2,2,2,1);Tin(1,2,1,1) Tin(1,2,2,2) Tin(1,2,1,2) Tin(1,2,2,1); Tin(2,1,1,1) Tin(2,1,2,2) Tin(2,1,1,2) Tin(2,1,2,1)];
T = [Spk(1,1) 0 Spk(1,2) 0; 0 Spk(2,2) 0 Spk(2,1); Spk(2,1) 0 Spk(2,2) 0; 0 Spk(1,2) 0 Spk(1,1)]; 
Small_sigma = Spk*F';

Sigma = [Small_sigma(1,1); Small_sigma(2,2); Small_sigma(2,1); Small_sigma(1,2)];

end
    

