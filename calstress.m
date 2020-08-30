function [stress] = calstress(eps,eta,x,y,disp_ele_left,Csmat, kron_delta, lamda,mu)
    F = zeros(3,3);
    [N,J,B] = shape(eps,eta,x,y);
    d_total = B*disp_ele_left;                
    F = [d_total(1,1) d_total(3,1) 0; d_total(4,1) d_total(2,1) 0; 0 0 0] + kron_delta;
    C = F'*F;
    E = 0.5*(C - kron_delta);
    [Dun,Tun,stress] = DTS(F,E,lamda,mu);
end