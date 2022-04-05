clc;
clear all;
% close all;
mu=.05;%learning rate form norm 2 factor
epsi=250;%l1 norm parameter 250
ro=.0007;%l1 norm parameter .00072
sample=1500;%input index
% burst=5;
for itr=1:200 %run times
    
    for j=1:5 %selecting different lms
        
     x=random('Normal',0,1,1,sample);
     noise=awgn(x,30)-x;
%      noise(300:300+burst-1)=80*rand(1,burst);
%      noise(1300:1300+burst-1)=80*rand(1,burst);
     sys=[zeros(1,4) 1 zeros(1,11)]; %system for input index 1:500
     sys_tap=zeros(1,length(sys)); %initialization
     model=zeros(1,length(sys));% model for the system
     model_tap=zeros(1,length(model));
     
    for i=1:500 %input index 1:500
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
%                 c(x2)=sign(model(x2)).*(1-19.*abs(model(x2)))./(1+abs(model(x2))).^21;%denominator power .5
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+.0001*c(x2); % new-sparse-lms
        end
        mdl4(i,:)=model;
        
    else 
                for x2=1:length(model)
                b=5;
                if(model(x2)>(-1/(b-1)) && model(x2)<(1/(b-1)))
                c(x2)=sign(model(x2)).*(1-4.*abs(model(x2)))./(1+abs(model(x2))).^5;%denominator power .5
%                     c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)-ro*c(x2); % new-sparse-lms
                end  
        mdl5(i,:)=model;
    end
    end    
    
    sys=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];  %system for input index 501:1000
     sys_tap=zeros(1,length(sys));
     model=zeros(1,length(sys));
     model_tap=zeros(1,length(model));
     
    for i=501:1000  %input index 501:1000
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
%                 c(x2)=sign(model(x2)).*(1-19.*abs(model(x2)))./(1+abs(model(x2))).^21;%denominator power .5
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+.0001*c(x2); % new-sparse-lms
        end
        mdl4(i,:)=model;
        
        
    else 
                for x2=1:length(model)
                b=5;
                if(model(x2)>(-1/(b-1)) && model(x2)<(1/(b-1)))
                c(x2)=sign(model(x2)).*(1-4.*abs(model(x2)))./(1+abs(model(x2))).^5;%denominator power .5
%                     c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)-ro*c(x2); % new-sparse-lms
                end
        
        mdl5(i,:)=model;
    end
    end

    sys=[1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];  %system for input index 1001:1500
     sys_tap=zeros(1,length(sys));
     model=zeros(1,length(sys));
     model_tap=zeros(1,length(model));
     
    for i=1001:1500 %input index 1001:1500
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
%                 c(x2)=sign(model(x2)).*(1-19.*abs(model(x2)))./(1+abs(model(x2))).^21;%denominator power .5
                    c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)+.0001*c(x2); % new-sparse-lms
        end
        mdl4(i,:)=model;
        
        
    else 
                for x2=1:length(model)
                b=5;
                if(model(x2)>(-1/(b-1)) && model(x2)<(1/(b-1)))
                c(x2)=sign(model(x2)).*(1-4.*abs(model(x2)))./(1+abs(model(x2))).^5;%denominator power .5
%                     c(x2)=((b)^2*(model(x2))-(b)*sign(model(x2)));
                else
                    c(x2)=0;
                end
                model(x2)=model(x2)+mu*(err(i))*model_tap(x2)-ro*c(x2); % new-sparse-lms
                end
        
        mdl5(i,:)=model;
    end
    end
  end %(j loop end)
     
for s=1:500  % MSD calculation at each input index
E1(itr,s)=mean(([zeros(1,4) 1 zeros(1,11)]-mdl1(s,:)).^2);
E2(itr,s)=mean(([zeros(1,4) 1 zeros(1,11)]-mdl2(s,:)).^2);
E3(itr,s)=mean(([zeros(1,4) 1 zeros(1,11)]-mdl3(s,:)).^2);
E4(itr,s)=mean(([zeros(1,4) 1 zeros(1,11)]-mdl4(s,:)).^2);
E5(itr,s)=mean(([zeros(1,4) 1 zeros(1,11)]-mdl5(s,:)).^2);
end
for s=501:1000;
E1(itr,s)=mean(([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-mdl1(s,:)).^2);
E2(itr,s)=mean(([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-mdl2(s,:)).^2);
E3(itr,s)=mean(([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-mdl3(s,:)).^2);
E4(itr,s)=mean(([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-mdl4(s,:)).^2);
E5(itr,s)=mean(([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]-mdl5(s,:)).^2);
end
for s=1001:1500;
E1(itr,s)=mean(([1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]-mdl1(s,:)).^2);
E2(itr,s)=mean(([1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]-mdl2(s,:)).^2);
E3(itr,s)=mean(([1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]-mdl3(s,:)).^2);
E4(itr,s)=mean(([1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]-mdl4(s,:)).^2);
E5(itr,s)=mean(([1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]-mdl5(s,:)).^2);
end
itr
end %(itr loop end)
figure;
plot(1:length(x),10*log10(mean(E1)),'--',1:length(x),10*log10(mean(E2)),'-',1:length(x),10*log10(mean(E3)),':',1:length(x),10*log10(mean(E4)),1:length(x),10*log10(mean(E5)),'Linewidth',2);
legend('lms','zalms','rzalms','l0-lms','new-sparse-lms');
xlabel('input index');
ylabel('MSD');
title('MSD plot ');


