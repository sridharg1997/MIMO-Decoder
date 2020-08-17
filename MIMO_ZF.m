function x_zf = MIMO_ZF(x,H,invtype)


[Nr,Nt] = size(H);
[~,N] = size(x);
temp = zeros(size(x));
for i =1:N
    y = zeros(Nr,Nt+1);
    y(:,1) = x(:,i);
    for k = 0:Nt-1
        e = zeros(Nt-k,1);
        e(Nt-k,1) = 1;
        %w = H(:,1:Nt-k)*((H(:,1:Nt-k)'*H(:,1:Nt-k))\e);
        w = H(:,1:Nt-k)*((H(:,1:Nt-k)'*H(:,1:Nt-k))\e);
        s_hat = w'*y(:,k+1);
        temp(Nr-k,i) = s_hat;
        y(:,k+2) = y(:,k+1) - H(:,Nt-k)*s_hat;
    end
end
x_zf = temp;