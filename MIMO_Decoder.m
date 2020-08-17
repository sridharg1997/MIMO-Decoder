%Tx-Rx setup 
N = 4;

Ncl = 6;
Nray = 8;
Nscatter = Nray*Ncl;


nsym = 2^12;
Rmod = 2;                         
M = 2^Rmod;
nbits = nsym*Rmod;


Niter = 100;

snr_param = -20:2:20;
Nsnr = numel(snr_param);
ber = zeros(Nsnr,2);

H = MIMO_Channel(N,N,Ncl,Nray);

for i = 1:Nsnr
    SNR = snr_param(i);                                % SNR in dB
    snr = db2pow(SNR);
    for j = 1:Niter

        bits = randi([0,1],nbits,1);
     

        sym = qammod(bits,M,'InputType','bit','UnitAveragePower',true);
        
        x = reshape(sym,[N,nsym/N]);

        Eb = 10*log10(mean(abs(sym).^2));

        N0_db = Eb - SNR;
        N0 = 10.^(0.1*N0_db);
        W =  sqrt(N0/2).*(randn(N,nsym/N) + 1i*randn(N,nsym/N));

        y = H*x + W;

        x_zf = MIMO_ZF(y,H);
        x_zf = reshape(x_zf,[nsym,1]);
        
        x_mmse = MIMO_MMSE(y,H,1/snr);
        x_mmse = reshape(x_mmse,[nsym,1]);
        %x_mmse =  MIMO_MMSE(y,H);

        bit_hat_zf = qamdemod(x_zf,M,'OutputType','bit');
        bit_hat_mmse = qamdemod(x_mmse,M,'OutputType','bit');

        ber(i,1) = ber(i,1) + mean(bit_hat_zf ~= bits);
        ber(i,2) = ber(i,2) + mean(bit_hat_mmse ~= bits);

    end
end

ber = ber/Niter;

semilogy(snr_param,ber(:,1),'b',snr_param,ber(:,2),'r');
grid on;
legend('ZF','MMSE');
xlabel('SNR in dB');
ylabel('BER');












