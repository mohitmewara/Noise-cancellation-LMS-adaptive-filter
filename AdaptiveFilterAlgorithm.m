input = importdata('project1.mat');
fs = input.fs;
primary = input.reference;
reference = input.primary;
primary_size = size(primary,2);
Epsilon = 0.0001;
 AllData = zeros(1,3);
 

for order = 30
    W_1 = zeros(order,1);  
    W_2 = zeros(order,1);
    
    performance_curve1 = zeros(46500,1);    
    performance_curve2 = zeros(18461,1);
    primary_wrt_filter = primary(1 , order:end);  %truncate primary
    reference_wrt_filter = zeros((primary_size - order),order);
    for update = (order) : primary_size                     %make reference_wrt_filter according to filter
        for update1=1:order
         reference_wrt_filter((update-order+1),update1) =  reference(update-update1+1);
        end  
    end
    
    disp(size(reference_wrt_filter,1));
    
%    for Nu = 0.001:(2 * 0.01):1
    for Nu = 0.05
%        for iterateReference = 1: size(reference_wrt_filter,1)
         for iterateReference = 1: size(reference_wrt_filter,1)
            MSE =0;
            Error = primary_wrt_filter(1, iterateReference) - (reference_wrt_filter(iterateReference,:) * W_2(:,1));
            X = reference_wrt_filter(iterateReference,:);
            Nu_by_Epsilon = Nu / (Epsilon + (X * X'));
            if iterateReference < 46501
                Error = primary_wrt_filter(1, iterateReference) - (reference_wrt_filter(iterateReference,:) * W_1(:,1));
                W_1 = W_1 + (Nu_by_Epsilon * (Error * X)');
                
                errorSquare = (primary_wrt_filter(1, 1:iterateReference)' - (reference_wrt_filter(1:iterateReference, :) * W_1(:,1))).^2;

                MSE = sum(errorSquare)/(iterateReference);
                performance_curve1(iterateReference,1) = MSE; 
                
            end  
            W_2 = W_2 + (Nu_by_Epsilon * (Error * X)');
            
            if iterateReference >= 46501
                errorSquare1 = (primary_wrt_filter(1, 46501:iterateReference)' - (reference_wrt_filter(46501:iterateReference, :) * W_2(:,1))).^2;

                MSE = sum(errorSquare1)/(iterateReference-46500);
                performance_curve2(iterateReference-46500,1) = MSE;
            end
                    
         end
%         MSE = sum(primary_wrt_filter - (reference_wrt_filter * W_1)')^2;
%         AllData(size(AllData,1)+1,1) = order;        
%         AllData(size(AllData,1),2) = Nu;
%         AllData(size(AllData,1),3) = MSE;
    end
    
Out = (primary_wrt_filter(1, 1:46500) - (reference_wrt_filter(1:46500,:) * W_1)');
Out1 = (primary_wrt_filter(1, 46501:end) - (reference_wrt_filter(46501:end,:) * W_2)');
Out3 = vertcat(Out', Out1');

SNR_parameter = mean(primary_wrt_filter.^2)/mean(Out3.^2);
SNR_After = 10 * log10(SNR_parameter);
figure;
plot(performance_curve1);
title('Learning Curve For Filter Order = 50 and Iteration < 46.5K');
xlabel('Iteration -->');
ylabel('MSE -->');
legend('Nu = 0.05');

figure;
plot(performance_curve2);
title('Learning Curve For Filter Order = 50 and Iteration > 46.5K');
xlabel('Iteration -->');
ylabel('MSE -->');
legend('Nu = 0.05');
% 
figure;
plot(Out3);
title('Error Signal After Applying NLMS For Filter Order = 50');
xlabel('Iteration -->');
ylabel('Error (Desired Output Signal) -->');
legend('Nu = 0.05');
%soundsc(Out3,fs);
end
 