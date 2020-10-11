clc;clear;
%Problem-8
%Generate random binary data of length-100000
l = 100000;%data length
data = randi([0 1], 1, l);
% 4-FSK
sf1 = [1,0,0,0]; % 00
sf2 = [0,1,0,0]; % 01
sf3 = [0,0,1,0]; % 10
sf4 = [0,0,0,1]; % 11
f=1;
for i = 1:2:length(data)%Encoding bits to symbols
    if data(i)==0 && data(i+1) == 0
        fsk(f,:)=sf1;f=f+1;
    elseif data(i)==0 && data(i+1) == 1
        fsk(f,:)=sf2;f=f+1;
    elseif data(i)==1 && data(i+1) == 0
        fsk(f,:)=sf3;f=f+1;
    else
        fsk(f,:)=sf4;f=f+1;
    end       
end

%Add AWGN noise to this signal. Vary the noise variance to have SNR range between [-5dB, 10 dB].
snr = -5:1:10;
for i=1:length(snr)
    rho = snr(i);
    j = (l/2)*(i-1)+1;
    n_f(j:j+((l/2)-1),1:4) = awgn(fsk,rho);%Noisy signal
end
%decoding symbols
df = zeros(length(n_f),4);
for j = 1:length(n_f)
    [m,p] = max(n_f(j,:));
    for i = 1:length(n_f(1,:))
        if(i==p)
            df(j,i)=1;
        else
            df(j,i)=0;
        end
    end    
end
%Symbol error rate
j=1;fsym_err = zeros(1,length(snr));
for s = 1:length(snr)
    fsym_err(s) = nnz(any(df(j:j+(l/2)-1,:)-fsk,2))/l;
    j = j+(l/2);
end
%Bit decoding from symbols
i=1;ber_fsk = zeros(1,l);
for b = 1:length(df)
    if(df(b,:)==[1 0 0 0])
        ber_fsk(i) = 0;
        ber_fsk(i+1)=0;
    elseif(df(b,:)==[0 1 0 0])
        ber_fsk(i) = 0;
        ber_fsk(i+1)=1;
    elseif(df(b,:)==[0 0 1 0])
        ber_fsk(i) = 1;
        ber_fsk(i+1)=0;
    else
        ber_fsk(i) = 1;
        ber_fsk(i+1)=1;
    end
    i = i+2;
end
%Bit error rate
j=1;fsk_ber = zeros(1,length(snr));
for s = 1:length(snr)
    fsk_ber(s) = nnz(ber_fsk(1,j:j+l-1)-data)/l;
    j = j+l;
end
%Plot
figure;
subplot(1,2,1);
plot(snr,fsym_err,'linewidth',2);
xlabel('SNR(in dB)');
ylabel('Symbol error rate(SER)');
title('SER vs SNR(in 4-FSK)');
grid;
set(gca,"linewidth",1,"fontsize",10);
subplot(1,2,2);
plot(snr,fsk_ber,'linewidth',2);
xlabel('SNR(in dB)');
ylabel('Bit error rate(BER)');
title('BER vs SNR(in 4-FSK)');
grid;
set(gca,"linewidth",1,"fontsize",10);