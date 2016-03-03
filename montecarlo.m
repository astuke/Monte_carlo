function mc()

    StepSize=linspace(0.1,1.0,10);
    NSteps=1000000;
    NThermalization=1000;
    
    O=zeros(length(StepSize),1);
    C1=zeros(length(StepSize),1);
    AcceptanceRatio=zeros(length(StepSize),1);
    ErrorEstimate=zeros(length(StepSize),1);
    
    for k=1:length(StepSize) %loop over step sizes
        
        [C1(k),AcceptanceRatio(k),O(k),ErrorEstimate(k)]=mc(StepSize(k),NSteps,NThermalization);
        
    end
    
    figure;
    plot(StepSize,C1);
    xlabel('Step size');
    ylabel('C(1)');
    
    figure;
    plot(StepSize,AcceptanceRatio);
    xlabel('Step size');
    ylabel('Acceptance ratio');
    
    figure;
    errorbar(StepSize,O,ErrorEstimate);
    xlabel('Step size');
    ylabel('O');
    

end

function [C1,AcceptanceRatio,ExpectationValue,ErrorEstimate]=mc(StepSize,NSteps,NThermalization)

    r=rand(6,1);
    
    for k=1:NThermalization
        
        [r,~]=update(r,StepSize);
        
    end
    
    OValues=zeros(NSteps,1);
    AcceptanceRatio=0;
    
    for k=1:NSteps
        
        [r,accepted]=update(r,StepSize);
        OValues(k)=operator(r); %operator value of current step 
        AcceptanceRatio=AcceptanceRatio+accepted;
        
    end
    
    AcceptanceRatio=AcceptanceRatio/NSteps;
    
    ExpectationValue=mean(OValues); %calculate expectation value of k operator values
    
    C=autocorr(OValues,1);
    C1=C(2);
    Nc=-1/log(C1)
    ErrorEstimate=std(OValues)/sqrt(NSteps/Nc); %estimate the error
    
    fprintf('StepSize=%f\n',StepSize);
    fprintf('<O>=%f +- %f\n',ExpectationValue,ErrorEstimate);
    fprintf('Acceptance ratio: %f\n',AcceptanceRatio);
    
end

function [r,accepted]=update(r,StepSize) %update function for importance sampling
    
    NElem=size(r);
    step=2*(rand(NElem(1),NElem(2))-0.5)*StepSize; % generates random number 
    rnew=r+step; %move electron position by random number -->update
    
    %accept update with probability
    if weight(r+step)/weight(r)>rand()
        r=r+step;
        accepted=1;
    else
        accepted=0;
    end

end

function O=operator(r)
    
    O=1/norm(r(1:3)-r(4:6)); %Coulomb operator
    %O=norm(r(1:3)); %|r_1|
    %O=norm(r(4:6)); %|r_2|
    %O=r(1); %x_1
    %O=r(4); %x_2
    
end

function g=weight(r)
    
    g=exp(-2*norm(r(1:3))-2*norm(r(4:6)))^2;
    
end
