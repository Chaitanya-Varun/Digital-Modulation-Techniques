clc;clear;
%Problem-8
%Generate random binary data of length-100000
l = 100000;%data length
data = randi([0 1], 1, l);
%BPSK
%Map the data to the BPSK signal constellation
Es = 1;%Energy of the signal
root_Es = sqrt(Es);
bpsk = (2*(data-0.5)*root_Es) ;%Mapping to constellation
%Add AWGN noise to this signal. Vary the noise variance to have SNR range between [-5dB, 10 dB].
snr = -5:1:10;
%wc = 1000;%carrier omega
for i=1:length(snr)
    rho = snr(i);
    n(i,:) = awgn(bpsk,rho);%Noisy signal
end
thresh_bpsk = 0;%Threshhold
for j = 1:length(snr)%Decoding
    for k = 1:length(bpsk)
        if(n(j,k)>thresh_bpsk)
            d(j,k) = 1;
        else
            d(j,k) = 0;
        end
    end
end
%Bit error rate
for m = 1:length(snr)
    ber_bpsk(m) = nnz(d(m,:)-data)/l;
end
%Plot
figure;
plot(snr,ber_bpsk,'linewidth',2);
xlabel('SNR(in dB)');
ylabel('Bit error rate(BER)');
title('BER vs SNR(in BPSK)');
grid;
set(gca,"linewidth",1,"fontsize",10);
