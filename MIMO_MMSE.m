function x_mmse = MIMO_MMSE(x,H,rho)

[Nr,Nt] = size(H);
[~,N] = size(x);
temp = zeros(size(x));
for i =1:N
    y = zeros(Nr,Nt+1);
    y(:,1) = x(:,i);
    for k = 0:Nr-1
        
        w = (H(:,1:Nr-k)*H(:,1:Nr-k)' + rho*eye(Nr))\H(:,Nt-k);
        s_hat = w'*y(:,k+1);
        temp(Nt-k,i) = s_hat;
        y(:,k+2) = y(:,k+1) - H(:,Nt-k)*s_hat;
    end
end
x_mmse = temp;