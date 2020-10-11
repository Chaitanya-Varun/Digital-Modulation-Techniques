%Question-2c
%Error function vs priory probability
p1 = 0:0.01:1;
sigma = sqrt(0.1); %assumed value
Ao = 1;%assumed value
for i = 1:101
    gamma = ((sigma*sigma)/(2*Ao))*log(p1(i)/(1-p1(i)));
    error(i) = p1(i)*(0.5*erfc((gamma+Ao)/(sigma*sqrt(2))))+ (1-p1(i))*(0.5*erfc((Ao-gamma)/(sigma*sqrt(2))));
end
    
plot(p1,error,'linewidth',2);
xlabel('Probability of occurence of symbol-1(p1)');
ylabel('Probability of bit error( Pr(error) )');
title('P(error) vs p1');
grid;
set(gca,"linewidth",1,"fontsize",10);