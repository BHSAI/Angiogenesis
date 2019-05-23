% Function to calculate the activation time, peah height and Deactivation time from all model output response curves

function [AT PH RI RP] = timescalculation(Y,tspan)
threshold = 0.5;% 50% of the peak value
Y = Y';
[F I] = max(Y);% find the value and location of the peak height of all model variables
PH = F;%Assign max values to PH 
AT = tspan(I);%Assign the location of peak height value as activation time for all model variables
Y50 = threshold.*F;% Calcuate 50% of the peak height values for all model variables

for i=1:length(F)

[R,C,V] = find(Y(:,i)>=Y50(i)); %Find the time location for 50% of peak height values
d = max(R);
DT(i) = tspan(d);
RI(i) = DT(i)-AT(i);% Resolution interval is the difference between the activation time and time at 50% of peak height values
end

RP = (Y(end,:)./F).*100;







