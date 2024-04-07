function [AIntra, fvalIntra] = GetATernarySixth(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';
X = ones(84, size(ROIpositions,1));

X(1,:) = ROIpositions(:,1);
X(2,:) = (ROIpositions(:,1).^2);
X(3,:) = (ROIpositions(:,1).^3);
X(4,:) = (ROIpositions(:,1).^4);
X(5,:) = (ROIpositions(:,1).^5);
X(6,:) = (ROIpositions(:,1).^6);

X(7,:) = ROIpositions(:,2);
X(8,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(9,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2);
X(10,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2);
X(11,:) = (ROIpositions(:,1).^4).*ROIpositions(:,2);
X(12,:) = (ROIpositions(:,1).^5).*ROIpositions(:,2);

X(13,:) = (ROIpositions(:,2).^2);
X(14,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2);
X(15,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2);
X(16,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,2).^2);
X(17,:) = (ROIpositions(:,1).^4).*(ROIpositions(:,2).^2);

X(18,:) = (ROIpositions(:,2).^3);
X(19,:) = ROIpositions(:,1).*(ROIpositions(:,2).^3);
X(20,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^3);
X(21,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,2).^3);

X(22,:) = (ROIpositions(:,2).^4);
X(23,:) = ROIpositions(:,1).*(ROIpositions(:,2).^4);
X(24,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^4);

X(25,:) = (ROIpositions(:,2).^5);
X(26,:) = ROIpositions(:,1).*(ROIpositions(:,2).^5);

X(27,:) = (ROIpositions(:,2).^6);

X(28,:) = ROIpositions(:,3);
X(29,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(30,:) = (ROIpositions(:,1).^2).*ROIpositions(:,3);
X(31,:) = (ROIpositions(:,1).^3).*ROIpositions(:,3);
X(32,:) = (ROIpositions(:,1).^4).*ROIpositions(:,3);
X(33,:) = (ROIpositions(:,1).^5).*ROIpositions(:,3);
X(34,:) = ROIpositions(:,2).*ROIpositions(:,3);
X(35,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);
X(36,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*ROIpositions(:,3);
X(37,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2).*ROIpositions(:,3);
X(38,:) = (ROIpositions(:,1).^4).*ROIpositions(:,2).*ROIpositions(:,3);
X(39,:) = (ROIpositions(:,2).^2).*ROIpositions(:,3);
X(40,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(41,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(42,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(43,:) = (ROIpositions(:,2).^3).*ROIpositions(:,3);
X(44,:) = ROIpositions(:,1).*(ROIpositions(:,2).^3).*ROIpositions(:,3);
X(45,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^3).*ROIpositions(:,3);
X(46,:) = (ROIpositions(:,2).^4).*ROIpositions(:,3);
X(47,:) = ROIpositions(:,1).*(ROIpositions(:,2).^4).*ROIpositions(:,3);
X(48,:) = (ROIpositions(:,2).^5).*ROIpositions(:,3);

X(49,:) = (ROIpositions(:,3).^2);
X(50,:) = ROIpositions(:,1).*(ROIpositions(:,3).^2);
X(51,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,3).^2);
X(52,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,3).^2);
X(53,:) = (ROIpositions(:,1).^4).*(ROIpositions(:,3).^2);
X(54,:) = ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(55,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(56,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(57,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(58,:) = (ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);
X(59,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);
X(60,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);
X(61,:) = (ROIpositions(:,2).^3).*(ROIpositions(:,3).^2);
X(62,:) = ROIpositions(:,1).*(ROIpositions(:,2).^3).*(ROIpositions(:,3).^2);
X(63,:) = (ROIpositions(:,2).^4).*(ROIpositions(:,3).^2);

X(64,:) = (ROIpositions(:,3).^3);
X(65,:) = ROIpositions(:,1).*(ROIpositions(:,3).^3);
X(66,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,3).^3);
X(67,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,3).^3);
X(68,:) = ROIpositions(:,2).*(ROIpositions(:,3).^3);
X(69,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^3);
X(70,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*(ROIpositions(:,3).^3);
X(71,:) = (ROIpositions(:,2).^2).*(ROIpositions(:,3).^3);
X(72,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*(ROIpositions(:,3).^3);
X(73,:) = (ROIpositions(:,2).^3).*(ROIpositions(:,3).^3);

X(74,:) = (ROIpositions(:,3).^4);
X(75,:) = ROIpositions(:,1).*(ROIpositions(:,3).^4);
X(76,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,3).^4);
X(77,:) = ROIpositions(:,2).*(ROIpositions(:,3).^4);
X(78,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^4);
X(79,:) = (ROIpositions(:,2).^2).*(ROIpositions(:,3).^4);

X(80,:) = (ROIpositions(:,3).^5);
X(81,:) = ROIpositions(:,1).*(ROIpositions(:,3).^5);
X(82,:) = ROIpositions(:,2).*(ROIpositions(:,3).^5);

X(83,:) = (ROIpositions(:,3).^6);

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];

% Subject to ...
Aeq = zeros(56,252);

Aeq(1,1) = 1; Aeq(1,91) = 1; Aeq(1,196) = 1; 
Aeq(2,2) = 2; Aeq(2,92) = 1; Aeq(2,197) = 1; 
Aeq(3,3) = 3; Aeq(3,93) = 1; Aeq(3,198) = 1; 
Aeq(4,4) = 4; Aeq(4,94) = 1; Aeq(4,199) = 1; 
Aeq(5,5) = 5; Aeq(5,95) = 1; Aeq(5,200) = 1; 
Aeq(6,6) = 6; Aeq(6,96) = 1; Aeq(6,201) = 1; 
Aeq(7,8) = 1; Aeq(7,97) = 2; Aeq(7,202) = 1;
Aeq(8,9) = 2; Aeq(8,98) = 2; Aeq(8,203) = 1; 
Aeq(9,10) = 3; Aeq(9,99) = 2; Aeq(9,204) = 1; 
Aeq(10,11) = 4; Aeq(10,100) = 2; Aeq(10,205) = 1; 
Aeq(11,12) = 5; Aeq(11,101) = 2; Aeq(11,206) = 1; 
Aeq(12,14) = 1; Aeq(12,102) = 3; Aeq(12,207) = 1; 
Aeq(13,15) = 2; Aeq(13,103) = 3; Aeq(13,208) = 1;
Aeq(14,16) = 3; Aeq(14,104) = 3; Aeq(14,209) = 1; 
Aeq(15,17) = 4; Aeq(15,105) = 3; Aeq(15,210) = 1; 
Aeq(16,19) = 1; Aeq(16,106) = 4; Aeq(16,211) = 1; 
Aeq(17,20) = 2; Aeq(17,107) = 4; Aeq(17,212) = 1; 
Aeq(18,21) = 3; Aeq(18,108) = 4; Aeq(18,213) = 1; 
Aeq(19,23) = 1; Aeq(19,109) = 5; Aeq(19,214) = 1; 
Aeq(20,24) = 2; Aeq(20,110) = 5; Aeq(20,215) = 1; 
Aeq(21,26) = 1; Aeq(21,111) = 6; Aeq(21,216) = 1; 
Aeq(22,29) = 1; Aeq(22,118) = 1; Aeq(22,217) = 2;
Aeq(23,30) = 2; Aeq(23,119) = 1; Aeq(23,218) = 2; 
Aeq(24,31) = 3; Aeq(24,120) = 1; Aeq(24,219) = 2; 
Aeq(25,32) = 4; Aeq(25,121) = 1; Aeq(25,220) = 2; 
Aeq(26,33) = 5; Aeq(26,122) = 1; Aeq(26,221) = 2;
Aeq(27,35) = 1; Aeq(27,123) = 2; Aeq(27,222) = 2; 
Aeq(28,36) = 2; Aeq(28,124) = 2; Aeq(28,223) = 2; 
Aeq(29,37) = 3; Aeq(29,125) = 2; Aeq(29,224) = 2; 
Aeq(30,38) = 4; Aeq(30,126) = 2; Aeq(30,225) = 2; 
Aeq(31,40) = 1; Aeq(31,127) = 3; Aeq(31,226) = 2; 
Aeq(32,41) = 2; Aeq(32,128) = 3; Aeq(32,227) = 2; 
Aeq(33,42) = 3; Aeq(33,129) = 3; Aeq(33,228) = 2; 
Aeq(34,44) = 1; Aeq(34,130) = 4; Aeq(34,229) = 2; 
Aeq(35,45) = 2; Aeq(35,131) = 4; Aeq(35,230) = 2; 
Aeq(36,47) = 1; Aeq(36,132) = 5; Aeq(36,231) = 2; 
Aeq(37,50) = 1; Aeq(37,138) = 1; Aeq(37,232) = 3; 
Aeq(38,51) = 2; Aeq(38,139) = 1; Aeq(38,233) = 3; 
Aeq(39,52) = 3; Aeq(39,140) = 1; Aeq(39,234) = 3; 
Aeq(40,53) = 4; Aeq(40,141) = 1; Aeq(40,235) = 3; 
Aeq(41,55) = 1; Aeq(41,142) = 2; Aeq(41,236) = 3; 
Aeq(42,56) = 2; Aeq(42,143) = 2; Aeq(42,237) = 3; 
Aeq(43,57) = 3; Aeq(43,144) = 2; Aeq(43,238) = 3; 
Aeq(44,59) = 1; Aeq(44,145) = 3; Aeq(44,239) = 3; 
Aeq(45,60) = 2; Aeq(45,146) = 3; Aeq(45,240) = 3; 
Aeq(46,62) = 1; Aeq(46,147) = 4; Aeq(46,241) = 3; 
Aeq(47,65) = 1; Aeq(47,152) = 1; Aeq(47,242) = 4; 
Aeq(48,66) = 2; Aeq(48,153) = 1; Aeq(48,243) = 4; 
Aeq(49,67) = 3; Aeq(49,154) = 1; Aeq(49,244) = 4; 
Aeq(50,69) = 1; Aeq(50,155) = 2; Aeq(50,245) = 4; 
Aeq(51,70) = 2; Aeq(51,156) = 2; Aeq(51,246) = 4; 
Aeq(52,72) = 1; Aeq(52,157) = 3; Aeq(52,247) = 4; 
Aeq(53,75) = 1; Aeq(53,161) = 1; Aeq(53,248) = 5; 
Aeq(54,76) = 2; Aeq(54,162) = 1; Aeq(54,249) = 5; 
Aeq(55,78) = 1; Aeq(55,163) = 2; Aeq(55,250) = 5; 
Aeq(56,81) = 1; Aeq(56,166) = 1; Aeq(56,251) = 6; 

Beq = zeros(56,1);

% Optimation
options = optimoptions('fmincon','Display','off');
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
AIntra = reshape(AIntra, 84, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20)-X(21,i)*x(21)-X(22,i)*x(22)-X(23,i)*x(23)-X(24,i)*x(24)-X(25,i)*x(25)-X(26,i)*x(26)-X(27,i)*x(27)-X(28,i)*x(28)-X(29,i)*x(29)-X(30,i)*x(30)-X(31,i)*x(31)-X(32,i)*x(32)-X(33,i)*x(33)-X(34,i)*x(34)-X(35,i)*x(35)-X(36,i)*x(36)-X(37,i)*x(37)-X(38,i)*x(38)-X(39,i)*x(39)-X(40,i)*x(40)-X(41,i)*x(41)-X(42,i)*x(42)-X(43,i)*x(43)-X(44,i)*x(44)-X(45,i)*x(45)-X(46,i)*x(46)-X(47,i)*x(47)-X(48,i)*x(48)-X(49,i)*x(49)-X(50,i)*x(50)-X(51,i)*x(51)-X(52,i)*x(52)-X(53,i)*x(53)-X(54,i)*x(54)-X(55,i)*x(55)-X(56,i)*x(56)-X(57,i)*x(57)-X(58,i)*x(58)-X(59,i)*x(59)-X(60,i)*x(60)-X(61,i)*x(61)-X(62,i)*x(62)-X(63,i)*x(63)-X(64,i)*x(64)-X(65,i)*x(65)-X(66,i)*x(66)-X(67,i)*x(67)-X(68,i)*x(68)-X(69,i)*x(69)-X(70,i)*x(70)-X(71,i)*x(71)-X(72,i)*x(72)-X(73,i)*x(73)-X(74,i)*x(74)-X(75,i)*x(75)-X(76,i)*x(76)-X(77,i)*x(77)-X(78,i)*x(78)-X(79,i)*x(79)-X(80,i)*x(80)-X(81,i)*x(81)-X(82,i)*x(82)-X(83,i)*x(83)-X(84,i)*x(84))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(84+1)-X(2,i)*x(84+2)-X(3,i)*x(84+3)-X(4,i)*x(84+4)-X(5,i)*x(84+5)-X(6,i)*x(84+6)-X(7,i)*x(84+7)-X(8,i)*x(84+8)-X(9,i)*x(84+9)-X(10,i)*x(84+10)-X(11,i)*x(84+11)-X(12,i)*x(84+12)-X(13,i)*x(84+13)-X(14,i)*x(84+14)-X(15,i)*x(84+15)-X(16,i)*x(84+16)-X(17,i)*x(84+17)-X(18,i)*x(84+18)-X(19,i)*x(84+19)-X(20,i)*x(84+20)-X(21,i)*x(84+21)-X(22,i)*x(84+22)-X(23,i)*x(84+23)-X(24,i)*x(84+24)-X(25,i)*x(84+25)-X(26,i)*x(84+26)-X(27,i)*x(84+27)-X(28,i)*x(84+28)-X(29,i)*x(84+29)-X(30,i)*x(84+30)-X(31,i)*x(84+31)-X(32,i)*x(84+32)-X(33,i)*x(84+33)-X(34,i)*x(84+34)-X(35,i)*x(84+35)-X(36,i)*x(84+36)-X(37,i)*x(84+37)-X(38,i)*x(84+38)-X(39,i)*x(84+39)-X(40,i)*x(84+40)-X(41,i)*x(84+41)-X(42,i)*x(84+42)-X(43,i)*x(84+43)-X(44,i)*x(84+44)-X(45,i)*x(84+45)-X(46,i)*x(84+46)-X(47,i)*x(84+47)-X(48,i)*x(84+48)-X(49,i)*x(84+49)-X(50,i)*x(84+50)-X(51,i)*x(84+51)-X(52,i)*x(84+52)-X(53,i)*x(84+53)-X(54,i)*x(84+54)-X(55,i)*x(84+55)-X(56,i)*x(84+56)-X(57,i)*x(84+57)-X(58,i)*x(84+58)-X(59,i)*x(84+59)-X(60,i)*x(84+60)-X(61,i)*x(84+61)-X(62,i)*x(84+62)-X(63,i)*x(84+63)-X(64,i)*x(84+64)-X(65,i)*x(84+65)-X(66,i)*x(84+66)-X(67,i)*x(84+67)-X(68,i)*x(84+68)-X(69,i)*x(84+69)-X(70,i)*x(84+70)-X(71,i)*x(84+71)-X(72,i)*x(84+72)-X(73,i)*x(84+73)-X(74,i)*x(84+74)-X(75,i)*x(84+75)-X(76,i)*x(84+76)-X(77,i)*x(84+77)-X(78,i)*x(84+78)-X(79,i)*x(84+79)-X(80,i)*x(84+80)-X(81,i)*x(84+81)-X(82,i)*x(84+82)-X(83,i)*x(84+83)-X(84,i)*x(84+84))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(168+1)-X(2,i)*x(168+2)-X(3,i)*x(168+3)-X(4,i)*x(168+4)-X(5,i)*x(168+5)-X(6,i)*x(168+6)-X(7,i)*x(168+7)-X(8,i)*x(168+8)-X(9,i)*x(168+9)-X(10,i)*x(168+10)-X(11,i)*x(168+11)-X(12,i)*x(168+12)-X(13,i)*x(168+13)-X(14,i)*x(168+14)-X(15,i)*x(168+15)-X(16,i)*x(168+16)-X(17,i)*x(168+17)-X(18,i)*x(168+18)-X(19,i)*x(168+19)-X(20,i)*x(168+20)-X(21,i)*x(168+21)-X(22,i)*x(168+22)-X(23,i)*x(168+23)-X(24,i)*x(168+24)-X(25,i)*x(168+25)-X(26,i)*x(168+26)-X(27,i)*x(168+27)-X(28,i)*x(168+28)-X(29,i)*x(168+29)-X(30,i)*x(168+30)-X(31,i)*x(168+31)-X(32,i)*x(168+32)-X(33,i)*x(168+33)-X(34,i)*x(168+34)-X(35,i)*x(168+35)-X(36,i)*x(168+36)-X(37,i)*x(168+37)-X(38,i)*x(168+38)-X(39,i)*x(168+39)-X(40,i)*x(168+40)-X(41,i)*x(168+41)-X(42,i)*x(168+42)-X(43,i)*x(168+43)-X(44,i)*x(168+44)-X(45,i)*x(168+45)-X(46,i)*x(168+46)-X(47,i)*x(168+47)-X(48,i)*x(168+48)-X(49,i)*x(168+49)-X(50,i)*x(168+50)-X(51,i)*x(168+51)-X(52,i)*x(168+52)-X(53,i)*x(168+53)-X(54,i)*x(168+54)-X(55,i)*x(168+55)-X(56,i)*x(168+56)-X(57,i)*x(168+57)-X(58,i)*x(168+58)-X(59,i)*x(168+59)-X(60,i)*x(168+60)-X(61,i)*x(168+61)-X(62,i)*x(168+62)-X(63,i)*x(168+63)-X(64,i)*x(168+64)-X(65,i)*x(168+65)-X(66,i)*x(168+66)-X(67,i)*x(168+67)-X(68,i)*x(168+68)-X(69,i)*x(168+69)-X(70,i)*x(168+70)-X(71,i)*x(168+71)-X(72,i)*x(168+72)-X(73,i)*x(168+73)-X(74,i)*x(168+74)-X(75,i)*x(168+75)-X(76,i)*x(168+76)-X(77,i)*x(168+77)-X(78,i)*x(168+78)-X(79,i)*x(168+79)-X(80,i)*x(168+80)-X(81,i)*x(168+81)-X(82,i)*x(168+82)-X(83,i)*x(168+83)-X(84,i)*x(168+84))*norm(Nvalues(i)))^2;
end
    
end

end

