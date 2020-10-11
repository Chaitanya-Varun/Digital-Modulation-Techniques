clc;clear;
%Problem-8
%Generate random binary data of length-100000
L = 100000;%data length
data = randi([0 1], 1, L);
%16-QAM
m = 16; K = sqrt(m);%M-ary
%Constellation
Q = zeros(m,2);
i=1;
for k = 1:K 
    x1 = 2*k-1-K;
    for l = 1:K
        x2 = 2*l-1-K;
        Q(i,:)= [ x1 x2 ];
        i=i+1;
    end
end
i=1;f=1;
qam = zeros(L/4,2);
for i = 1:4:L-K+1
    if data(i)==0 && data(i+1)==1 && data(i+2)==0 && data(i+3)==0
        qam(f,:) = Q(1,:);
    elseif data(i)==1 && data(i+1)==1 && data(i+2)==0 && data(i+3)==0
        qam(f,:) = Q(2,:);
    elseif data(i)==1 && data(i+1)==0 && data(i+2)==0 && data(i+3)==0
        qam(f,:) = Q(3,:);
    elseif data(i)==0 && data(i+1)==0 && data(i+2)==0 && data(i+3)==0
        qam(f,:) = Q(4,:); 
    elseif data(i)==0 && data(i+1)==1 && data(i+2)==0 && data(i+3)==1
        qam(f,:) = Q(5,:);
    elseif data(i)==1 && data(i+1)==1 && data(i+2)==0 && data(i+3)==1
        qam(f,:) = Q(6,:);
    elseif data(i)==1 && data(i+1)==0 && data(i+2)==0 && data(i+3)==1
        qam(f,:) = Q(7,:);
    elseif data(i)==0 && data(i+1)==0 && data(i+2)==0 && data(i+3)==1
        qam(f,:) = Q(8,:);
    elseif data(i)==0 && data(i+1)==1 && data(i+2)==1 && data(i+3)==1
        qam(f,:) = Q(9,:);
    elseif data(i)==1 && data(i+1)==1 && data(i+2)==1 && data(i+3)==1
        qam(f,:) = Q(10,:);
    elseif data(i)==1 && data(i+1)==0 && data(i+2)==1 && data(i+3)==1
        qam(f,:) = Q(11,:);
    elseif data(i)==0 && data(i+1)==0 && data(i+2)==1 && data(i+3)==1
        qam(f,:) = Q(12,:);
    elseif data(i)==0 && data(i+1)==1 && data(i+2)==1 && data(i+3)==0
        qam(f,:) = Q(13,:);
    elseif data(i)==1 && data(i+1)==1 && data(i+2)==1 && data(i+3)==0
        qam(f,:) = Q(14,:);
    elseif data(i)==1 && data(i+1)==0 && data(i+2)==1 && data(i+3)==0
        qam(f,:) = Q(15,:);
    else
        qam(f,:) = Q(16,:);
    end
    f = f+1;
end
%Adding noise to the signal
snr = -5:1:10;qam_n = zeros(4*L,2);
for i=1:length(snr)
    rho = snr(i);
    j = (L/4)*(i-1)+1;
    qam_n(j:(j-1+(L/4)),:) = awgn(qam,rho);%noisy signal
end
%Decoding of symbol
i=1;
qam_d = zeros(L*4,2);
for i=1:L*4
    for d=1:2
        if(qam_n(i,d)>2)
            qam_d(i,d) = 3;
        elseif(qam_n(i,d)>=0 && qam_n(i,d)<=2)
            qam_d(i,d) = 1;
        elseif(qam_n(i,d)<=0 && qam_n(i,d)>=-2)
            qam_d(i,d) = -1;
        else
            qam_d(i,d) = -3;
        end
    end
end
%Symbol error rate
j=1;qsym_err = zeros(1,length(snr));
for s = 1:length(snr)
    qsym_err(s) = nnz(any(qam_d(j:j+(L/4)-1,:)-qam,2))/l;
    j = j+(L/4);
end

%Converting symbols to bits
i=1;p=1;bqam= zeros(L,1);
for i = 1:length(qam_d)
    if(qam_d(i,:)==Q(1,:))
        bqam(p)=0;bqam(p+1)=1;bqam(p+2)=0;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(2,:))
        bqam(p)=1;bqam(p+1)=1;bqam(p+2)=0;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(3,:))
        bqam(p)=1;bqam(p+1)=0;bqam(p+2)=0;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(4,:))
        bqam(p)=0;bqam(p+1)=0;bqam(p+2)=0;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(5,:))
        bqam(p)=0;bqam(p+1)=1;bqam(p+2)=0;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(6,:))
        bqam(p)=1;bqam(p+1)=1;bqam(p+2)=0;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(7,:))
        bqam(p)=1;bqam(p+1)=0;bqam(p+2)=0;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(8,:))
        bqam(p)=0;bqam(p+1)=0;bqam(p+2)=0;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(9,:))
        bqam(p)=0;bqam(p+1)=1;bqam(p+2)=1;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(10,:))
        bqam(p)=1;bqam(p+1)=1;bqam(p+2)=1;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(11,:))
        bqam(p)=1;bqam(p+1)=0;bqam(p+2)=1;bqam(p+3)=1;    
    elseif(qam_d(i,:)==Q(12,:))
        bqam(p)=0;bqam(p+1)=0;bqam(p+2)=1;bqam(p+3)=1;
    elseif(qam_d(i,:)==Q(13,:))
        bqam(p)=0;bqam(p+1)=1;bqam(p+2)=1;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(14,:))
        bqam(p)=1;bqam(p+1)=1;bqam(p+2)=1;bqam(p+3)=0;
    elseif(qam_d(i,:)==Q(15,:))
        bqam(p)=1;bqam(p+1)=0;bqam(p+2)=1;bqam(p+3)=0;
    else
         bqam(p)=0;bqam(p+1)=0;bqam(p+2)=1;bqam(p+3)=0;
    end
    p=p+4;
end
%bit error rate
j=1;qam_ber = zeros(1,length(snr));
bqam = reshape(bqam,[1 16*L]);
for s = 1:length(snr)
    qam_ber(s) = nnz(bqam(1,j:j+L-1)-data)/l;
    j = j+L;
end
%Plot
figure;
subplot(1,2,1);
plot(snr,qsym_err,'linewidth',2);
xlabel('SNR(in dB)');
ylabel('Symbol error rate(SER)');
title('SER vs SNR(in 16-QAM)');
grid;
set(gca,"linewidth",1,"fontsize",10);
subplot(1,2,2);
plot(snr,qam_ber,'linewidth',2);
xlabel('SNR(in dB)');
ylabel('Bit error rate(BER)');
title('BER vs SNR(in 16-QAM)');
grid;
set(gca,"linewidth",1,"fontsize",10);











