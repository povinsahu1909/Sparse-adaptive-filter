clc;
clear all;
close all;
format long

function [h]=rir(fs, mic, n, r, rm, src);
%RIR   Room Impulse Response.
%   [h] = RIR(FS, MIC, N, R, RM, SRC) performs a room impulse
%         response calculation by means of the mirror image method.
%
%      FS  = sample rate.
%      MIC = row vector giving the x,y,z coordinates of
%            the microphone.  
%      N   = The program will account for (2*N+1)^3 virtual sources 
%      R   = reflection coefficient for the walls, in general -1<R<1.
%      RM  = row vector giving the dimensions of the room.  
%      SRC = row vector giving the x,y,z coordinates of 
%            the sound source.
%
%   EXAMPLE:
% %
%      fs=44100;
%      mic=[19 18 1.6];
%      n=12;
%      r=0.3;
%      rm=[20 19 21];
%      src=[5 2 1];
%      rir(fs, mic, n, r, rm, src);
%
%   NOTES:
%
%   1) All distances are in meters.
%   2) The output is scaled such that the largest value of the 
%      absolute value of the output vector is equal to one.
%   3) To implement this filter, you will need to do a fast 
%      convolution.  The program FCONV.m will do this. It can be 
%      found on the Mathworks File Exchange at
%      www.mathworks.com/matlabcentral/fileexchange/.  It can also 
%      be found at http://www.sgm-audio.com/research/rir/fconv.m
%   4) A paper has been written on this model.  It is available at:
%      http://www.sgm-audio.com/research/rir/rir.html
%      
%
%Version 3.4.2
%Copyright © 2003 Stephen G. McGovern

%Some of the following comments are references to equations the my paper.

nn=-n:1:n;                            % Index for the sequence
rms=nn+0.5-0.5*(-1).^nn;              % Part of equations 2,3,& 4
srcs=(-1).^(nn);                      % part of equations 2,3,& 4
xi=srcs*src(1)+rms*rm(1)-mic(1);      % Equation 2 
yj=srcs*src(2)+rms*rm(2)-mic(2);      % Equation 3 
zk=srcs*src(3)+rms*rm(3)-mic(3);      % Equation 4 

[i,j,k]=meshgrid(xi,yj,zk);           % convert vectors to 3D matrices
d=sqrt(i.^2+j.^2+k.^2);               % Equation 5
time=round(fs*d/343)+1;               % Similar to Equation 6
              
[e,f,g]=meshgrid(nn, nn, nn);         % convert vectors to 3D matrices
c=r.^(abs(e)+abs(f)+abs(g));          % Equation 9
e=c./d;                               % Equivalent to Equation 10

h=full(sparse(time(:),1,e(:)));       % Equivalent to equation 11
h=h/max(abs(h));                       % Scale output
end

% DONE take u as 0.0006 and p as 0.00004
fs=500;
mic=[19 18 1.6];
n=8;
r=0.3;
rm=[20 19 21];
src=[7 8 1];
sys2=(rir(fs, mic, n, r, rm, src))';

% figure;
% plot(sys3);
count3=0;
for i=1:length(sys2)
    if sys2(i)<0.00005
        count3=count3+1;
    end
end
figure;
plot(sys2)


%parameters for robust error
sigma_square=2;

%parameters to calculate ANR(Average signal to noise ratio)
A_e=0;
A_d=0;
lambda=0.85;

for n=1:2
    n
    sys=sys2;
%input sequence
x=rand(1,25000)-0.5;
% x=filter([1 2 1],[1 -1.45 0.57],x);

noise=awgn(x,40)-x;
noise(8000:8100)=100*rand(1,101);
noise(16000:16100)=100*rand(1,101);
dim=(length(sys)*2);
dim1=length(sys);
y=[zeros(1,floor(dim1/2)),x];
x_cap=zeros(1,dim1);
tap=zeros(1,length(sys));
x_cap1=zeros(1,dim1);
tap1=zeros(1,length(sys));
 
%sparse parameters
dim_PNLMS=(length(sys));
w_sparse=rand(1,dim_PNLMS);
sys_tap_sparse=zeros(1,length(sys));
sys_tap_model_sparse=zeros(1,dim_PNLMS);
u_sparse=0.08;

%sparse IPNMLS parameters
w_IPsparse=rand(1,dim_PNLMS);
sys_tap_IPsparse=zeros(1,length(sys));
sys_tap_model_IPsparse=zeros(1,dim_PNLMS);

%parameters for PNLMS
sys_tap_PNLMS=zeros(1,length(sys));
% u_PNLMS=0.2;
% u_PNLMS=2.4;
u_PNLMS=2;
sys_tap_model_PNLMS=zeros(1,dim_PNLMS);
w_PNLMS=rand(1,dim_PNLMS);
alpha1 =0;
alpha2=0;
% %lms robust
% w_robust=rand(1,dim);
% sys_tap_robust=zeros(1,length(sys));
% sys_tap_model_robust=zeros(1,dim);
% u_robust=0.001;

%parameters for LMS
sys_tap=zeros(1,length(sys));
sys_tap_model=zeros(1,dim1);
w=rand(1,dim1);
% u=0.1;
u=9e-02;


%parameters to calculate ANR(Average signal to noise ratio)
A_e=0;
A_d=0;
A_e_robust=0;
A_d_robust=0;
A_e_PNLMS=0;
A_d_PNLMS=0;
A_e_sparse=0;
A_d_sparse=0;
lambda=0.85;
    for i=1:length(x)
        
        sys_tap_model=[x(i),sys_tap_model(1:end-1)];
        ModelValue(i)=sys_tap_model*w';
        sys_tap=[ModelValue(i),sys_tap(1:end-1)];
        output(i)=sys_tap*sys';
        e(i)=y(i)-output(i)+noise(i);
        
%         %robust
%         sys_tap_model_robust=[x(i),sys_tap_model_robust(1:end-1)];
%         ModelValue_robust(i)=sys_tap_model_robust*w_robust';
%         sys_tap_robust=[ModelValue_robust(i),sys_tap_robust(1:end-1)];
%         output_robust(i)=sys_tap_robust*sys';
%         e_robust(i)=y(i)-output_robust(i)+noise(i);
        
        sys_tap_model_PNLMS=[x(i),sys_tap_model_PNLMS(1:end-1)];
        ModelValue_PNLMS(i)=sys_tap_model_PNLMS*w_PNLMS';
        sys_tap_PNLMS=[ModelValue_PNLMS(i),sys_tap_PNLMS(1:end-1)];
        output_PNLMS(i)=sys_tap_PNLMS*sys';
        e_PNLMS(i)=y(i)-output_PNLMS(i)+noise(i);
        
        sys_tap_model_sparse=[x(i),sys_tap_model_sparse(1:end-1)];
        ModelValue_sparse(i)=sys_tap_model_sparse*w_sparse';
        sys_tap_sparse=[ModelValue_sparse(i),sys_tap_sparse(1:end-1)];
        output_sparse(i)=sys_tap_sparse*sys';
        e_sparse(i)=y(i)-output_sparse(i)+noise(i);
        
        %to calculate x_cap
        tap=[x(i),tap(1:end-1)];
        sys_opt(i)=tap*sys';
        x_cap = [sys_opt(i),x_cap(1:end-1)];
        
        %to calculate x_cap
        tap1=[x(i),tap1(1:end-1)];
        sys_opt(i)=tap1*sys';
        x_cap1 = [sys_opt(i),x_cap1(1:end-1)];
        
        %lms
        w=w+2*u*e(i)*x_cap1/(x_cap1*x_cap1'+0.001);
        
%         %lms robust
%         w_robust=w_robust+2*u_robust*(e_robust(i)/((e_robust(i).^2)+2*sigma_square))*x_cap;
            
        %robust IPNLMS
        norm_h=sum(w_PNLMS);
        delta = (1-alpha2)*var(x)/(2*length(w_PNLMS));
        delta_cap = (1-alpha2)*var(x_cap)/(2*length(w_PNLMS));
        for j=1:length(w_PNLMS)
            k(j)= ((1-alpha2)*norm_h/length(w_PNLMS))+(1+alpha2)*abs(w_PNLMS(j));
        end
% %         a(i,:)=k;
        K = diag(k);
        if i<300
            w_PNLMS=w_PNLMS+2*(u_PNLMS*K*x_cap1'*e_PNLMS(i))'/((e_PNLMS(i).^2+2.5)*(x_cap1*K*x_cap1'+7));
        else
            w_PNLMS=w_PNLMS+2*(u_PNLMS*K*x_cap1'*e_PNLMS(i))'/((e_PNLMS(i).^2+2.5)*(x_cap1*K*x_cap1'+7));
        end
            %         if i<101
%             w_PNLMS=w_PNLMS+2*(u_PNLMS*K*x_cap'*e_PNLMS(i))'/((e_PNLMS(i).^2+2.5)*(x_cap*K*x_cap'+4));
%         else    
%             w_PNLMS=w_PNLMS+2*(u_PNLMS*K*x_cap'*e_PNLMS(i))'/((e_PNLMS(i).^2+0.1*std(e_PNLMS(i-100:i)))*(x_cap*K*x_cap'+4));
%         end
        
        %IPNLMS
        norm_h=sum(w_sparse);
        delta = (1-alpha1)*var(x)/(2*length(w_sparse));
        delta_cap = (1-alpha1)*var(x_cap1)/(2*length(w_sparse));
        for j=1:length(w_sparse)
            t(j)= ((1-alpha1)*norm_h/length(w_sparse))+(1+alpha1)*abs(w_sparse(j));
        end

        T = diag(t);
        if i<300
            w_sparse=w_sparse+(2*u_sparse*T*x_cap1'*e_sparse(i)/(x_cap1*T*x_cap1'+var(x)))';   
        else
            w_sparse=w_sparse+(2*u_sparse*T*x_cap1'*e_sparse(i)/(x_cap1*T*x_cap1'+var(x)))';
        end
            A_e = lambda*A_e + (1-lambda)*abs(e(i));
            A_d = lambda*A_d + (1-lambda)*abs(y(i));
            ANR(i) = 20*log10(A_e/A_d);
        
%             A_e_robust = lambda*A_e_robust + (1-lambda)*abs(e_robust(i));
%             A_d_robust = lambda*A_d_robust + (1-lambda)*abs(y(i));
%             ANR_robust(i) = 20*log10(A_e_robust/A_d_robust);

            A_e_PNLMS = lambda*A_e_PNLMS + (1-lambda)*abs(e_PNLMS(i));
            A_d_PNLMS = lambda*A_d_PNLMS + (1-lambda)*abs(y(i));
            ANR_PNLMS(i) = 20*log10(A_e_PNLMS/A_d_PNLMS);
        
            A_e_sparse = lambda*A_e_sparse + (1-lambda)*abs(e_sparse(i));
            A_d_sparse = lambda*A_d_sparse + (1-lambda)*abs(y(i));
            ANR_sparse(i) = 20*log10(A_e_sparse/A_d_sparse);
        
            
            if i==17000
               store(1,:)=w;
               store(2,:)=w_PNLMS;
               store(3,:)=w_sparse;
            end
    end    
    err(n,:)=e.^2;
    err_IPNLMS(n,:)=e_PNLMS.^2;
    err_sparse(n,:)=e_sparse.^2;
%     err_robust(n,:)=e_robust.^2;
    anr_lms(n,:)=ANR;
    anr_PNLMS(n,:)=ANR_PNLMS;
    anr_lms_sparse(n,:)=ANR_sparse;
%     anr_robust(n,:)=ANR_robust;
end    
%     w
% e(end).^2
% plot(ANR)
% figure
% plot(store)

figure
plot(y);
hold on
plot(output_PNLMS,'r');
hold on
plot(output,'g');

figure
plot(10*log10(mean(err)),'b');
hold on
plot(10*log10(mean(err_IPNLMS)),'r');
hold on
plot(10*log10(mean(err_sparse)),'y');
% hold on
% plot(10*log10(mean(err_robust)),'y');

% figure
% plot(10*log10(e.^2),'b');
% hold on
% plot(10*log10(e_PNLMS.^2),'r');
% hold on
% plot(10*log10(e_sparse.^2),'g');

figure
plot(mean(anr_lms),'b');
hold on
plot(mean(anr_PNLMS),'r');
hold on
plot(mean(anr_lms_sparse),'y');
% hold on
% plot(10*log10(mean(anr_robust)),'y');

% figure;
% freqz(sys)
% figure;
% freqz(store(1,:),'r')
% figure;
% freqz(store(2,:),'g')
% figure;
% freqz(store(3,:),'y')
% figure;
% a=freqz(sys)+freqz(store(1,:));
% plot(:,1)
% hold on
% plot(10*log10(e_IPsparse.^2),'c');
% plot(y);
% hold on
% plot(output,'r')
s=[1:512]./512;
 figure;
 plot(s,20*log10(abs(freqz(store(1,:))))+20*log10(abs(freqz(sys))))
 hold on
 plot(s,20*log10(abs(freqz(store(2,:))))+20*log10(abs(freqz(sys))))
 hold on
 plot(s,20*log10(abs(freqz(store(3,:))))+20*log10(abs(freqz(sys))),'r')