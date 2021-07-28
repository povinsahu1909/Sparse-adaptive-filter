clc;
clear all;
%%
inp_len=10000;
for i=1:100
    i
x=rand(1,inp_len)-0.5;% input
h=[-0.5 0.1 1 0.1 -0.5]; % unknown system
d=filter(h,1,x);% unknown system output
d=awgn(d,40);
num_of_coeff=10;
x_tap=zeros(1,num_of_coeff);
w=zeros(1,num_of_coeff);
mu=0.02;
for n=1:inp_len
    x_tap=[x(n) x_tap(1:end-1)];
    y(n)=w*x_tap';
    e(n)=d(n)-y(n);
    w=w+mu*e(n)*x_tap;
end

Err(i,:)=e.^2;
end
%%
figure
plot(10*log10(mean(Err)))

    
    
    
    
    