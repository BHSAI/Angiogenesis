function [ATf_KO PHf_KO DTf_KO RPf_KO Y_output] = knockoutcalculation(Parametersets_impaired_secondfilter_TGFKO,tmax,tspan,yinit,delaytime,delaytime2) 
  
 % parfor_progress(iter)
  parpool('local',4)
parfor i = 1:size(Parametersets_impaired_secondfilter_TGFKO,1)
%for i = 1:size(Parametersets_impaired_secondfilter_TGFKO,1)
    [AT1 PH1 DT1 RP1 Y_main1]= Curvecharacteristic(tmax,tspan,yinit,Parametersets_impaired_secondfilter_TGFKO(i,:),delaytime,delaytime2);
    ATf_KO(i,:) = AT1;
    PHf_KO(i,:) = PH1;
    DTf_KO(i,:) = DT1;
    RPf_KO(i,:) = RP1;
    Y_output{i} = Y_main1;   
end
delete(gcp('nocreate'));