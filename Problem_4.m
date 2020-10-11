% Problem-4 
n=6;%number of bits per sample
L=64;%levels
f_s = 2* 4000;%Sampling frequency
T_s = 1/f_s;%Time for each sample
T_b = T_s/6;%Time for each bit
%Instruction-1 : first generate a random PCM sequence of length 100
% As we got 6-bits per sample, the numbers we can use is in range [0,63]
%rng(0,'twister');
ri = 63*rand(1,100);
ro = floor(ri);%Quantizing
qe_max = 1;
q_np = qe_max^2/3*(4^n);
r = ro;
%r = randi([0 63],1,100);
E = norm(r)/100;%Average enegy of signal
sqnr = E/q_np;
sqnr_db = 10*log(sqnr)/log(10);
b = de2bi(r');
bn = b(:).'; % Now we are having 600 bits 
%Applying 4-PAM to this bits
%Assuming the levels for '00'-> -3;'01'-> -1;'10'-> 1;'11'-> 3;
j=1;%Iteration purposes
for i = 1:2:length(bn)-1
    if bn(i)==0 && bn(i+1)==0
        pam(j)=-3;j = j+1;
    elseif bn(i)==0 && bn(i+1)==1
        pam(j)=-1;j = j+1;
    elseif bn(i)==1 && bn(i+1)==0
        pam(j)=1;j = j+1;
    elseif bn(i)==1 && bn(i+1)==1
        pam(j)=3;j = j+1;
    end
end
%Now we have the pam signals
%According to instruction 2 we need to sampled at 4times the designed baseband sample frequency
k=0;
for x = 1:length(pam)
    s_pam(k+1) = pam(x);
    s_pam(k+2) = pam(x);
    s_pam(k+3) = pam(x);
    s_pam(k+4) = pam(x);
    k = k+4;
end
T_sym = 2*T_b;
%Now we have the transmitted waveform
%From instruction-3, we have to add AWGN noise so that the SNR is 6 dB at
%the input of the matched filter, 
E_s = (5/T_sym)*T_sym;
%The general method 
%snr = 6;
%sigma = sqrt(E_s*10^(-snr/10));%Here sigma = 1.12;
%mu = 0;%the mean os gaussian white noise is 0;
%n = sigma*randn(size(s_pam))+mu;
%out = s_pam+n;
out = awgn(s_pam,6,10*log(5)/log(10));%All the above functions are done with AWGN
%plot(x.*T_sym,out)
%Instruction-4: Implement the matched filter and sampler at the channel sample frequency
%Implementing matched filter
m = [ 1 1 1 1];
z1 = filter(m,1,out);%Function which does match filter
%The convolution method
%z= [0 0 0 0 0 0 0];
%for i=1:4:length(out)-3;
 %z(1,i:i+6) = [z(i:i+2) 0 0 0 0]+conv(out(1,i:i+3),[1 1 1 1]);%convolution of pulses with m(t)
 %end
 %Normalising
 z1 = z1./sqrt(20);
 
%Now sampling at channel frequency
for i = 1:length(z1)
    if mod(i,4) == 0
        s = i/4;
        z_s(s)=z1(i);
    end
end
%Here z_s contains all the values of signal sampled at k*T_sym 
%Instruction-5: Make a decision at the output of the Matched filter.
dec_pcm = [];%decode PCM
dec_pam = [];%decode PAM
for e = 1:length(z_s)
    m = min([(z_s(e)+3)*(z_s(e)+3),(z_s(e)+1)*(z_s(e)+1),(z_s(e)-1)*(z_s(e)-1),(z_s(e)-3)*(z_s(e)-3)]);
    if m==(z_s(e)+3)*(z_s(e)+3)
        dec_pcm=cat(2,dec_pcm,[0 0]);
        dec_pam=cat(2,dec_pam,(-3));
    elseif m==(z_s(e)+1)*(z_s(e)+1)
        dec_pcm=cat(2,dec_pcm,[0 1]);
        dec_pam=cat(2,dec_pam,(-1));
    elseif m==(z_s(e)-1)*(z_s(e)-1)
        dec_pcm=cat(2,dec_pcm,[1 0]);
        dec_pam=cat(2,dec_pam,(1));
    else 
        dec_pcm=cat(2,dec_pcm,[1 1]);
        dec_pam=cat(2,dec_pam,(3));
    end     
end
dec_b = reshape(dec_pcm,[100,6]);% decoded binary symbols
%Instruction 6 : Compare the transmitted and decoded PAM symbols and find the symbol error rate.
%Transmitted PAM symbol is given by pam[], decoded PAM symbols is dec_pam[]
sym_err = nnz(dec_pam-pam);% Number errors in the decoded symbols
sym_err_rate = (sym_err/length(pam))*100;%symbol error rate
%Instruction 7 : find the PCM symbol error rate and bit error rate. What is the SQNR?
bit_err = nnz(dec_pcm-bn);%Number of errors in decoded PCM bits
bit_err_rate = (bit_err/length(dec_pcm))*100;%bit error rate
err_dec_sam = 100-sum( all( (dec_b - b) == 0, 2 ) ); % Number of wrongly decoded samples
disp("Number errors in the decoded symbols - "+sym_err+" Symbol error rate :" +sym_err_rate+"%");
disp("Number of errors in decoded PCM bits - "+bit_err+" Bit error rate :"+bit_err_rate+"%");
disp("Number of wrongly decoded samples - "+err_dec_sam);
disp("SQNR - "+sqnr+" SQNR in dB :"+sqnr_db+" dB");

%Plotting(Please open the figure in full screen)
x_1=1:1:600;x_1 = x_1.*T_b;
x_2 = 1:1:300;x_2=x_2.*4;
x = 1:1:1200;
subplot(4,1,1);
plot(x_1,bn,'linewidth',1.25);
xlim([ 0 0.006]);
title('PCM Sequence');
xlabel('Time(accordance with T_b)');
ylabel('Amplitude');
set(gca,"linewidth",1,"fontsize",10);
subplot(4,1,2);
plot(x,s_pam,'linewidth',1.25);
title('4-PAM Signal');
xlabel('Time scaled as T_{b}');
ylabel('Amplitude');
set(gca,"linewidth",1,"fontsize",10);
subplot(4,1,3);
plot(x,out,'linewidth',1.25);
title('4-PAM Signal+AWGN');
xlabel('Time scaled as T_{b}');
ylabel('Amplitude');
set(gca,"linewidth",1,"fontsize",10);
subplot(4,1,4);
plot(x,z1,'linewidth',1.5);
title('Matched Filter output,Sampled output');
xlabel('Time scaled as T_{b}');
ylabel('Amplitude');
hold on;
plot(x_2,z_s,'ro');
legend("Matched Filter output","Sampled output");
hold off;
set(gca,"linewidth",1,"fontsize",10);














