function [N,J, B] = shape(eps,eta,x,y)
            N = (1/4)*[(1-eps)*(1-eta),0,(1+eps)*(1-eta),0,(1+eps)*(1+eta),0,(1-eps)*(1+eta),0;0,(1-eps)*(1-eta),0,(1+eps)*(1-eta),0,(1+eps)*(1+eta),0,(1-eps)*(1+eta)];
            dx_deps = (1/4)*(-(1-eta)*x(1)+(1-eta)*x(2)+(1+eta)*x(3)-(1+eta)*x(4));
            dx_deta = (1/4)*(-(1-eps)*x(1)-(1+eps)*x(2)+(1+eps)*x(3)+(1-eps)*x(4));
            dy_deps = (1/4)*(-(1-eta)*y(1)+(1-eta)*y(2)+(1+eta)*y(3)-(1+eta)*y(4));
            dy_deta = (1/4)*(-(1-eps)*y(1)-(1+eps)*y(2)+(1+eps)*y(3)+(1-eps)*y(4));
            dN_deps = (1/4)*[-(1-eta),(1-eta),(1+eta),-(1+eta)];
            dN_deta = (1/4)*[-(1-eps),-(1+eps),(1+eps),(1-eps)];
            J = [dx_deps,dy_deps;dx_deta,dy_deta];
            
            Bnew = double(J\[dN_deps; dN_deta]); 
            B1 = [Bnew(1,1),0;0,Bnew(2,1);Bnew(2,1),0; 0,Bnew(1,1)];
            B2 = [Bnew(1,2),0;0,Bnew(2,2);Bnew(2,2),0;0,Bnew(1,2)];
            B3 = [Bnew(1,3),0;0,Bnew(2,3);Bnew(2,3),0;0,Bnew(1,3)];
            B4 = [Bnew(1,4),0;0,Bnew(2,4);Bnew(2,4),0;0,Bnew(1,4)];
            B = [B1,B2,B3,B4];
end
