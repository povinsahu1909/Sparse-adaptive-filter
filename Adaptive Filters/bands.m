clc;
clear all;
% close all;
mu=.05; %learning rate form norm 2 factor

sample=500; %input index
% x=audioread('Timit.wav');
% sample=length(x);

epsi=250; %l1 norm parameter
ro=.0007; %l1 norm parameter
rand('state',1);
for itr=1:200%run times
    
    for j=1:10 %selecting different lms

        x=random('Normal',0,1,1,sample);
        noise=awgn(x,30)-x;
%         noise=random('Normal',0,10^(-1.5),1,sample);
        sys=[zeros(1,4) 1 zeros(1,11)];
    % sys=[zeros(1,4) .9 -.9 .8 -.8 .7 -.7 .6 -.6 .5 -.5 .4 -.4 .3 -.3 .2 -.2 .1 -.1 zeros( 1,178) .5 -.5 .4 -.4 .3 -.3 .2 -.2 .1 -.1 zeros(1,48)]; %system for input index 1:5000
%     sys=[zeros(1,60) -0.0001 -0.0002 -0.002 -0.001 0.15 0.12  -0.14 -0.13 -0.12 0.05 0.065 0.045 0.03 0.012 0 -0.01 -0.03 -0.05 0.025 0.025 -0.02 -0.025 0.043 0.02 -0.05 -0.045 0.02 0.03 0.02 -0.01 0 0.025 -0.002 0.01 0 zeros(1,25)];
     sys_tap=zeros(1,length(sys)); %initialization
     model=zeros(1,length(sys));% model for the system
     model_tap=zeros(1,length(model));
     
    for i=1:sample %input index 1:500
    sys_tap=[x(i) sys_tap(1:end-1)];
    sys_out(i)=sys_tap*sys'+noise(i);
    model_tap=[x(i)  model_tap(1:end-1)];
    model_out(i)=model_tap*model';
    err(i)=sys_out(i)-model_out(i);
    if(j==1)
       model=model+mu*(err(i))*model_tap; %lms
       mdl1(i,:)=model;
    elseif(j==2)
        model=model+mu*(err(i))*model_tap-ro*sign(model); %za-lms
        mdl2(i,:)=model;
    elseif(j==3)
        model=model+mu*err(i)*model_tap-ro*sign(model)./(1+200*abs(model)); %RZA-lms
        mdl3(i,:)=model;
    elseif(j==4)
        for x2=1:length(model)
                b=4;
                if(model(x2)>(-1/(b)) && model(x2)<(1/(b)))
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+.0001*c(x2); % new-sparse-lms
        end
        mdl4(i,:)=model;
        
        
    elseif(j==5) 
        for x2=1:length(model)
                b=2;
                if(model(x2)>(-1/(b-1)) && model(x2)<(1/(b-1)))
                    c(x2)=sign(model(x2)).*(1-4.*abs(model(x2)))./(1+abs(model(x2))).^5;%denominator power .5
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)-ro*c(x2); % new-sparse-lms
        end
        mdl5(i,:)=model;

%     Two band ZALMS
    elseif(j==6)
        for x2=1:length(model)
                b=2;
                if (model(x2)>(-1/(2*b)) && model(x2)<(1/(2*b)))
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                elseif(model(x2)>(-1/(b)) && model(x2)<(1/(b)))
                    c(x2)=((b/2)^2*(model(x2))-(b/2)*sign(model(x2)));
                else 
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+ro*c(x2); % new-sparse-lms
        end
        mdl6(i,:)=model;
            
    % Two band RZALMS 
    elseif(j==7)
        for x2=1:length(model)
                b=3;
                if (model(x2)>(-1/(2*(b-1))) && model(x2)<(1/(2*(b-1))))
                    c(x2)=sign(model(x2)).*(1-(b-1).*abs(model(x2)))./(1+abs(model(x2))).^b;
                elseif(model(x2)>(-1/(b-1)) && model(x2)<(1/(b-1)))
                    c(x2)=sign(model(x2)).*(1-(2*(b-1)).*abs(model(x2)))./(1+abs(model(x2))).^(2*b);%denominator power .5
                else 
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)-ro*c(x2); % new-sparse-lms
        end
        mdl7(i,:)=model;

        %     Three band ZALMS
    elseif(j==8)
        for x2=1:length(model)
                b=3;
                if (model(x2)>(-1/(3*b)) && model(x2)<(1/(3*b)))
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                elseif (model(x2)>(-1/(2*b)) && model(x2)<(1/(2*b)))
                    c(x2)=((b/2)^2*(model(x2))-(b/2)*sign(model(x2)));
                elseif(model(x2)>(-1/(b)) && model(x2)<(1/(b)))
                    c(x2)=((b/4)^2*(model(x2))-(b/4)*sign(model(x2)));
                else 
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+ro*c(x2); % new-sparse-lms
        end
        mdl8(i,:)=model;
    

    %     Four band ZALMS
    elseif(j==9)
        for x2=1:length(model)
                b=4;
                if (model(x2)>(-1/(4*b)) && model(x2)<(1/(4*b)))
                        c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));    
                elseif (model(x2)>(-1/(3*b)) && model(x2)<(1/(3*b)))
                    c(x2)=((b/2)^2*(model(x2))-(b/2)*sign(model(x2)));
                elseif (model(x2)>(-1/(2*b)) && model(x2)<(1/(2*b)))
                    c(x2)=((b/4)^2*(model(x2))-(b/4)*sign(model(x2)));
                elseif(model(x2)>(-1/(b)) && model(x2)<(1/(b)))
                    c(x2)=((b/8)^2*(model(x2))-(b/8)*sign(model(x2)));
                else 
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+ro*c(x2); % new-sparse-lms
        end
        mdl9(i,:)=model;

    %     Five band ZALMS
    else
        for x2=1:length(model)
                b=5;
                if (model(x2)>(-1/(5*b)) && model(x2)<(1/(5*b)))
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2))); 
                elseif (model(x2)>(-1/(4*b)) && model(x2)<(1/(4*b)))
                    c(x2)=((b/2)^2*(model(x2))-(b/2)*sign(model(x2)));    
                elseif (model(x2)>(-1/(3*b)) && model(x2)<(1/(3*b)))
                    c(x2)=((b/4)^2*(model(x2))-(b/4)*sign(model(x2)));
                elseif (model(x2)>(-1/(2*b)) && model(x2)<(1/(2*b)))
                    c(x2)=((b/8)^2*(model(x2))-(b/8)*sign(model(x2)));
                elseif(model(x2)>(-1/(b)) && model(x2)<(1/(b)))
                    c(x2)=((b/16)^2*(model(x2))-(b/16)*sign(model(x2)));
                else 
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+ro*c(x2); % new-sparse-lms
        end
        mdl10(i,:)=model;

    end
    end
    error(j,:)=err.^2;
  end %(j loop end)
    
  err1(itr,:)=error(1,:);
  err2(itr,:)=error(2,:);
  err3(itr,:)=error(3,:);
  err4(itr,:)=error(4,:);
  err5(itr,:)=error(5,:);
  err6(itr,:)=error(6,:);
  err7(itr,:)=error(7,:);
  err8(itr,:)=error(8,:);
for s=1:sample  % MSD calculation at each input index
E1(itr,s)=mean((sys-mdl1(s,:)).^2);
E2(itr,s)=mean((sys-mdl2(s,:)).^2);
E3(itr,s)=mean((sys-mdl3(s,:)).^2);
E4(itr,s)=mean((sys-mdl4(s,:)).^2);
E5(itr,s)=mean((sys-mdl5(s,:)).^2);
E6(itr,s)=mean((sys-mdl6(s,:)).^2);
E7(itr,s)=mean((sys-mdl7(s,:)).^2);
E8(itr,s)=mean((sys-mdl8(s,:)).^2);
E9(itr,s)=mean((sys-mdl9(s,:)).^2);
E10(itr,s)=mean((sys-mdl10(s,:)).^2);
end
itr
end %(itr loop end)
% E1=(err1);
% E2=(err2);
% E3=(err3);
% E4=(err4);
figure;
plot(1:length(x),10*log10(mean(E1)),1:length(x),10*log10(mean(E2)),1:length(x),10*log10(mean(E3)),1:length(x),10*log10(mean(E4)),1:length(x),10*log10(mean(E5)),1:length(x),10*log10(mean(E6)), 1:length(x),10*log10(mean(E7)), 1:length(x),10*log10(mean(E8)), 1:length(x),10*log10(mean(E9)), 1:length(x),10*log10(mean(E10)));
legend('LMS','ZA-LMS','RZA-LMS','BB-ZALMS','BB-RZALMS', 'MBB-ZALMS', 'MBB-RZALMS', 'band3', 'band4', 'band5');
xlabel('Iterations');
ylabel('MSD (dB)');

figure;
plot(1:length(x),10*log10(mean(E4)),1:length(x),10*log10(mean(E5)),1:length(x),10*log10(mean(E6)), 1:length(x),10*log10(mean(E7)), 1:length(x),10*log10(mean(E8)), 1:length(x),10*log10(mean(E9)), 1:length(x),10*log10(mean(E10)));
legend('BB-ZALMS','BB-RZALMS', 'band2-ZALMS', 'band2-RZALMS', 'band3-ZALMS', 'band4-ZALMS', 'band5-ZALMS');
xlabel('Iterations');
ylabel('MSD (dB)');

% figure;
% plot(sys);
% hold on;
% plot(model,'r');
% e=10*log10(mean(E1));
