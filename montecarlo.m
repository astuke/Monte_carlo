function mc()
 NSteps=10000;
 r=rand(6,1);
 O=zeros(NSteps,1);
 acc=0;
 NLags=1;
 C1=zeros(0,1);
 
 for i=0.1:0.1:1   %this is the the loop over the step size
     
    for k=1:NSteps   % loop over all steps
    [r,acc]=update(r,acc,i);
    O(k)=operator(r);
    end
    
  C=autocorr(O,NLags); %autocorrelation function 
  C1=[C1 C(2)];
  %disp('Accpetance = ') 
  acc=acc/NSteps % calculate and output acceptance --> should be 50%
 end
 
 C1

 
 ExpectedO=mean(O) % expectation value of operator
 %hist(O,100); % error histogram
 %plot(C);
 sigma=std(O);
 Nc=(1+C(2))/(1-C(2));
 error=sigma/sqrt(NSteps/Nc)   % estimated error of result
end

function [r,acc]=update(r,acc,i)
 StepSize=i;
 delta=2*(rand(6,1)-0.5)*StepSize; %uniformly distributed random numbers from -1 to 1 multiplied with StepSize
    if rand()<weight(r+delta)/weight(r)  
    r=r+delta; %if accepted: add random value
    acc=acc+1;
    end  
end

function w=weight(r)
 w=exp(-2*norm(r(1:3))-2*norm(r(4:6)))^2;
end

function val=operator(r)
 val=1/norm(r(1:3)-r(4:6)); %value of Coulomb potential energy
 %val=norm(r(4:6)); %value of vector r1
 %val=r(2);   %value of single coordinate
 end
