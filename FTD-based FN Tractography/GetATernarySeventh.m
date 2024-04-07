function [AIntra, fvalIntra] = GetATernarySeventh(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';
X = ones(120, size(ROIpositions,1));

X(1,:) = ROIpositions(:,1);
X(2,:) = ROIpositions(:,1).^2;
X(3,:) = ROIpositions(:,1).^3;
X(4,:) = ROIpositions(:,1).^4;
X(5,:) = ROIpositions(:,1).^5;
X(6,:) = ROIpositions(:,1).^6;
X(7,:) = ROIpositions(:,1).^7;
X(8,:) = ROIpositions(:,2);
X(9,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(10,:) = ROIpositions(:,1).^2.*ROIpositions(:,2);
X(11,:) = ROIpositions(:,1).^3.*ROIpositions(:,2);
X(12,:) = ROIpositions(:,1).^4.*ROIpositions(:,2);
X(13,:) = ROIpositions(:,1).^5.*ROIpositions(:,2);
X(14,:) = ROIpositions(:,1).^6.*ROIpositions(:,2);
X(15,:) = ROIpositions(:,2).^2;
X(16,:) = ROIpositions(:,1).*ROIpositions(:,2).^2;
X(17,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2;
X(18,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2;
X(19,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^2;
X(20,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).^2;
X(21,:) = ROIpositions(:,2).^3;
X(22,:) = ROIpositions(:,1).*ROIpositions(:,2).^3;
X(23,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3;
X(24,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^3;
X(25,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^3;
X(26,:) = ROIpositions(:,2).^4;
X(27,:) = ROIpositions(:,1).*ROIpositions(:,2).^4;
X(28,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^4;
X(29,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^4;
X(30,:) = ROIpositions(:,2).^5;
X(31,:) = ROIpositions(:,1).*ROIpositions(:,2).^5;
X(32,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^5;
X(33,:) = ROIpositions(:,2).^6;
X(34,:) = ROIpositions(:,1).*ROIpositions(:,2).^6;
X(35,:) = ROIpositions(:,2).^7;
X(36,:) = ROIpositions(:,3);
X(37,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(38,:) = ROIpositions(:,1).^2.*ROIpositions(:,3);
X(39,:) = ROIpositions(:,1).^3.*ROIpositions(:,3);
X(40,:) = ROIpositions(:,1).^4.*ROIpositions(:,3);
X(41,:) = ROIpositions(:,1).^5.*ROIpositions(:,3);
X(42,:) = ROIpositions(:,1).^6.*ROIpositions(:,3);
X(43,:) = ROIpositions(:,2).*ROIpositions(:,3);
X(44,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);
X(45,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3);
X(46,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3);
X(47,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).*ROIpositions(:,3);
X(48,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).*ROIpositions(:,3);
X(49,:) = ROIpositions(:,2).^2.*ROIpositions(:,3);
X(50,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(51,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(52,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(53,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(54,:) = ROIpositions(:,2).^3.*ROIpositions(:,3);
X(55,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(56,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(57,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(58,:) = ROIpositions(:,2).^4.*ROIpositions(:,3);
X(59,:) = ROIpositions(:,1).*ROIpositions(:,2).^4.*ROIpositions(:,3);
X(60,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^4.*ROIpositions(:,3);
X(61,:) = ROIpositions(:,2).^5.*ROIpositions(:,3);
X(62,:) = ROIpositions(:,1).*ROIpositions(:,2).^5.*ROIpositions(:,3);
X(63,:) = ROIpositions(:,2).^6.*ROIpositions(:,3);
X(64,:) = ROIpositions(:,3).^2;
X(65,:) = ROIpositions(:,1).*ROIpositions(:,3).^2;
X(66,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^2;
X(67,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^2;
X(68,:) = ROIpositions(:,1).^4.*ROIpositions(:,3).^2;
X(69,:) = ROIpositions(:,1).^5.*ROIpositions(:,3).^2;
X(70,:) = ROIpositions(:,2).*ROIpositions(:,3).^2;
X(71,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(72,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(73,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(74,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(75,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(76,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(77,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(78,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(79,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(80,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(81,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(82,:) = ROIpositions(:,2).^4.*ROIpositions(:,3).^2;
X(83,:) = ROIpositions(:,1).*ROIpositions(:,2).^4.*ROIpositions(:,3).^2;
X(84,:) = ROIpositions(:,2).^5.*ROIpositions(:,3).^2;
X(85,:) = ROIpositions(:,3).^3;
X(86,:) = ROIpositions(:,1).*ROIpositions(:,3).^3;
X(87,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^3;
X(88,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^3;
X(89,:) = ROIpositions(:,1).^4.*ROIpositions(:,3).^3;
X(90,:) = ROIpositions(:,2).*ROIpositions(:,3).^3;
X(91,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(92,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(93,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(94,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(95,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(96,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(97,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^3;
X(98,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3).^3;
X(99,:) = ROIpositions(:,2).^4.*ROIpositions(:,3).^3;
X(100,:) = ROIpositions(:,3).^4;
X(101,:) = ROIpositions(:,1).*ROIpositions(:,3).^4;
X(102,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^4;
X(103,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^4;
X(104,:) = ROIpositions(:,2).*ROIpositions(:,3).^4;
X(105,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^4;
X(106,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^4;
X(107,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^4;
X(108,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^4;
X(109,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^4;
X(110,:) = ROIpositions(:,3).^5;
X(111,:) = ROIpositions(:,1).*ROIpositions(:,3).^5;
X(112,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^5;
X(113,:) = ROIpositions(:,2).*ROIpositions(:,3).^5;
X(114,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^5;
X(115,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^5;
X(116,:) = ROIpositions(:,3).^6;
X(117,:) = ROIpositions(:,1).*ROIpositions(:,3).^6;
X(118,:) = ROIpositions(:,2).*ROIpositions(:,3).^6;
X(119,:) = ROIpositions(:,3).^7;

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];

% Subject to ...
Aeq = zeros(84,360);

Aeq(1,1) = 1; Aeq(1,128) = 1; Aeq(1,276) = 1; 
Aeq(2,2) = 2; Aeq(2,129) = 1; Aeq(2,277) = 1; 
Aeq(3,3) = 3; Aeq(3,130) = 1; Aeq(3,278) = 1; 
Aeq(4,4) = 4; Aeq(4,131) = 1; Aeq(4,279) = 1; 
Aeq(5,5) = 5; Aeq(5,132) = 1; Aeq(5,280) = 1; 
Aeq(6,6) = 6; Aeq(6,133) = 1; Aeq(6,281) = 1; 
Aeq(7,7) = 7; Aeq(7,134) = 1; Aeq(7,282) = 1;
Aeq(8,9) = 1; Aeq(8,135) = 2; Aeq(8,283) = 1; 
Aeq(9,10) = 2; Aeq(9,136) = 2; Aeq(9,284) = 1; 
Aeq(10,11) = 3; Aeq(10,137) = 2; Aeq(10,285) = 1; 
Aeq(11,12) = 4; Aeq(11,138) = 2; Aeq(11,286) = 1; 
Aeq(12,13) = 5; Aeq(12,139) = 2; Aeq(12,287) = 1; 
Aeq(13,14) = 6; Aeq(13,140) = 2; Aeq(13,288) = 1;
Aeq(14,16) = 1; Aeq(14,141) = 3; Aeq(14,289) = 1; 
Aeq(15,17) = 2; Aeq(15,142) = 3; Aeq(15,290) = 1; 
Aeq(16,18) = 3; Aeq(16,143) = 3; Aeq(16,291) = 1; 
Aeq(17,19) = 4; Aeq(17,144) = 3; Aeq(17,292) = 1; 
Aeq(18,20) = 5; Aeq(18,145) = 3; Aeq(18,293) = 1; 
Aeq(19,22) = 1; Aeq(19,146) = 4; Aeq(19,294) = 1; 
Aeq(20,23) = 2; Aeq(20,147) = 4; Aeq(20,295) = 1; 
Aeq(21,24) = 3; Aeq(21,148) = 4; Aeq(21,296) = 1; 
Aeq(22,25) = 4; Aeq(22,149) = 4; Aeq(22,297) = 1;
Aeq(23,27) = 1; Aeq(23,150) = 5; Aeq(23,298) = 1; 
Aeq(24,28) = 2; Aeq(24,151) = 5; Aeq(24,299) = 1; 
Aeq(25,29) = 3; Aeq(25,152) = 5; Aeq(25,300) = 1; 
Aeq(26,31) = 1; Aeq(26,153) = 6; Aeq(26,301) = 1;
Aeq(27,32) = 2; Aeq(27,154) = 6; Aeq(27,302) = 1; 
Aeq(28,34) = 1; Aeq(28,155) = 7; Aeq(28,303) = 1; 
Aeq(29,37) = 1; Aeq(29,163) = 1; Aeq(29,304) = 2; 
Aeq(30,38) = 2; Aeq(30,164) = 1; Aeq(30,305) = 2; 
Aeq(31,39) = 3; Aeq(31,165) = 1; Aeq(31,306) = 2; 
Aeq(32,40) = 4; Aeq(32,166) = 1; Aeq(32,307) = 2; 
Aeq(33,41) = 5; Aeq(33,167) = 1; Aeq(33,308) = 2; 
Aeq(34,42) = 6; Aeq(34,168) = 1; Aeq(34,309) = 2; 
Aeq(35,44) = 1; Aeq(35,169) = 2; Aeq(35,310) = 2; 
Aeq(36,45) = 2; Aeq(36,170) = 2; Aeq(36,311) = 2; 
Aeq(37,46) = 3; Aeq(37,171) = 2; Aeq(37,312) = 2; 
Aeq(38,47) = 4; Aeq(38,172) = 2; Aeq(38,313) = 2; 
Aeq(39,48) = 5; Aeq(39,173) = 2; Aeq(39,314) = 2; 
Aeq(40,50) = 1; Aeq(40,174) = 3; Aeq(40,315) = 2; 
Aeq(41,51) = 2; Aeq(41,175) = 3; Aeq(41,316) = 2; 
Aeq(42,52) = 3; Aeq(42,176) = 3; Aeq(42,317) = 2; 
Aeq(43,53) = 4; Aeq(43,177) = 3; Aeq(43,318) = 2; 
Aeq(44,55) = 1; Aeq(44,178) = 4; Aeq(44,319) = 2; 
Aeq(45,56) = 2; Aeq(45,179) = 4; Aeq(45,320) = 2; 
Aeq(46,57) = 3; Aeq(46,180) = 4; Aeq(46,321) = 2; 
Aeq(47,59) = 1; Aeq(47,181) = 5; Aeq(47,322) = 2; 
Aeq(48,60) = 2; Aeq(48,182) = 5; Aeq(48,323) = 2; 
Aeq(49,62) = 1; Aeq(49,183) = 6; Aeq(49,324) = 2; 
Aeq(50,65) = 1; Aeq(50,190) = 1; Aeq(50,325) = 3; 
Aeq(51,66) = 2; Aeq(51,191) = 1; Aeq(51,326) = 3; 
Aeq(52,67) = 3; Aeq(52,192) = 1; Aeq(52,327) = 3; 
Aeq(53,68) = 4; Aeq(53,193) = 1; Aeq(53,328) = 3; 
Aeq(54,69) = 5; Aeq(54,194) = 1; Aeq(54,329) = 3; 
Aeq(55,71) = 1; Aeq(55,195) = 2; Aeq(55,330) = 3; 
Aeq(56,72) = 2; Aeq(56,196) = 2; Aeq(56,331) = 3; 
Aeq(57,73) = 3; Aeq(57,197) = 2; Aeq(57,332) = 3; 
Aeq(58,74) = 4; Aeq(58,198) = 2; Aeq(58,333) = 3; 
Aeq(59,76) = 1; Aeq(59,199) = 3; Aeq(59,334) = 3; 
Aeq(60,77) = 2; Aeq(60,200) = 3; Aeq(60,335) = 3; 
Aeq(61,78) = 3; Aeq(61,201) = 3; Aeq(61,336) = 3; 
Aeq(62,80) = 1; Aeq(62,202) = 4; Aeq(62,337) = 3; 
Aeq(63,81) = 2; Aeq(63,203) = 4; Aeq(63,338) = 3; 
Aeq(64,83) = 1; Aeq(64,204) = 5; Aeq(64,339) = 3; 
Aeq(65,86) = 1; Aeq(65,210) = 1; Aeq(65,340) = 4; 
Aeq(66,87) = 2; Aeq(66,211) = 1; Aeq(66,341) = 4; 
Aeq(67,88) = 3; Aeq(67,212) = 1; Aeq(67,342) = 4; 
Aeq(68,89) = 4; Aeq(68,213) = 1; Aeq(68,343) = 4; 
Aeq(69,91) = 1; Aeq(69,214) = 2; Aeq(69,344) = 4; 
Aeq(70,92) = 2; Aeq(70,215) = 2; Aeq(70,345) = 4; 
Aeq(71,93) = 3; Aeq(71,216) = 2; Aeq(71,346) = 4; 
Aeq(72,95) = 1; Aeq(72,217) = 3; Aeq(72,347) = 4; 
Aeq(73,96) = 2; Aeq(73,218) = 3; Aeq(73,348) = 4; 
Aeq(74,98) = 1; Aeq(74,219) = 4; Aeq(74,349) = 4; 
Aeq(75,101) = 1; Aeq(75,224) = 1; Aeq(75,350) = 5; 
Aeq(76,102) = 2; Aeq(76,225) = 1; Aeq(76,351) = 5; 
Aeq(77,103) = 3; Aeq(77,226) = 1; Aeq(77,352) = 5; 
Aeq(78,105) = 1; Aeq(78,227) = 2; Aeq(78,353) = 5; 
Aeq(79,106) = 2; Aeq(79,228) = 2; Aeq(79,354) = 5; 
Aeq(80,108) = 1; Aeq(80,229) = 3; Aeq(80,355) = 5; 
Aeq(81,111) = 1; Aeq(81,233) = 1; Aeq(81,356) = 6; 
Aeq(82,112) = 2; Aeq(82,234) = 1; Aeq(82,357) = 6; 
Aeq(83,114) = 1; Aeq(83,235) = 2; Aeq(83,358) = 6; 
Aeq(84,117) = 1; Aeq(84,238) = 1; Aeq(84,359) = 7; 

Beq = zeros(84,1);

% Optimation
options = optimoptions('fmincon','Display','off');
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
AIntra = reshape(AIntra, 120, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20)-X(21,i)*x(21)-X(22,i)*x(22)-X(23,i)*x(23)-X(24,i)*x(24)-X(25,i)*x(25)-X(26,i)*x(26)-X(27,i)*x(27)-X(28,i)*x(28)-X(29,i)*x(29)-X(30,i)*x(30)-X(31,i)*x(31)-X(32,i)*x(32)-X(33,i)*x(33)-X(34,i)*x(34)-X(35,i)*x(35)-X(36,i)*x(36)-X(37,i)*x(37)-X(38,i)*x(38)-X(39,i)*x(39)-X(40,i)*x(40)-X(41,i)*x(41)-X(42,i)*x(42)-X(43,i)*x(43)-X(44,i)*x(44)-X(45,i)*x(45)-X(46,i)*x(46)-X(47,i)*x(47)-X(48,i)*x(48)-X(49,i)*x(49)-X(50,i)*x(50)-X(51,i)*x(51)-X(52,i)*x(52)-X(53,i)*x(53)-X(54,i)*x(54)-X(55,i)*x(55)-X(56,i)*x(56)-X(57,i)*x(57)-X(58,i)*x(58)-X(59,i)*x(59)-X(60,i)*x(60)-X(61,i)*x(61)-X(62,i)*x(62)-X(63,i)*x(63)-X(64,i)*x(64)-X(65,i)*x(65)-X(66,i)*x(66)-X(67,i)*x(67)-X(68,i)*x(68)-X(69,i)*x(69)-X(70,i)*x(70)-X(71,i)*x(71)-X(72,i)*x(72)-X(73,i)*x(73)-X(74,i)*x(74)-X(75,i)*x(75)-X(76,i)*x(76)-X(77,i)*x(77)-X(78,i)*x(78)-X(79,i)*x(79)-X(80,i)*x(80)-X(81,i)*x(81)-X(82,i)*x(82)-X(83,i)*x(83)-X(84,i)*x(84)-X(85,i)*x(85)-X(86,i)*x(86)-X(87,i)*x(87)-X(88,i)*x(88)-X(89,i)*x(89)-X(90,i)*x(90)-X(91,i)*x(91)-X(92,i)*x(92)-X(93,i)*x(93)-X(94,i)*x(94)-X(95,i)*x(95)-X(96,i)*x(96)-X(97,i)*x(97)-X(98,i)*x(98)-X(99,i)*x(99)-X(100,i)*x(100)-X(101,i)*x(101)-X(102,i)*x(102)-X(103,i)*x(103)-X(104,i)*x(104)-X(105,i)*x(105)-X(106,i)*x(106)-X(107,i)*x(107)-X(108,i)*x(108)-X(109,i)*x(109)-X(110,i)*x(110)-X(111,i)*x(111)-X(112,i)*x(112)-X(113,i)*x(113)-X(114,i)*x(114)-X(115,i)*x(115)-X(116,i)*x(116)-X(117,i)*x(117)-X(118,i)*x(118)-X(119,i)*x(119)-X(120,i)*x(120))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(121)-X(2,i)*x(122)-X(3,i)*x(123)-X(4,i)*x(124)-X(5,i)*x(125)-X(6,i)*x(126)-X(7,i)*x(127)-X(8,i)*x(128)-X(9,i)*x(129)-X(10,i)*x(130)-X(11,i)*x(131)-X(12,i)*x(132)-X(13,i)*x(133)-X(14,i)*x(134)-X(15,i)*x(135)-X(16,i)*x(136)-X(17,i)*x(137)-X(18,i)*x(138)-X(19,i)*x(139)-X(20,i)*x(140)-X(21,i)*x(141)-X(22,i)*x(142)-X(23,i)*x(143)-X(24,i)*x(144)-X(25,i)*x(145)-X(26,i)*x(146)-X(27,i)*x(147)-X(28,i)*x(148)-X(29,i)*x(149)-X(30,i)*x(150)-X(31,i)*x(151)-X(32,i)*x(152)-X(33,i)*x(153)-X(34,i)*x(154)-X(35,i)*x(155)-X(36,i)*x(156)-X(37,i)*x(157)-X(38,i)*x(158)-X(39,i)*x(159)-X(40,i)*x(160)-X(41,i)*x(161)-X(42,i)*x(162)-X(43,i)*x(163)-X(44,i)*x(164)-X(45,i)*x(165)-X(46,i)*x(166)-X(47,i)*x(167)-X(48,i)*x(168)-X(49,i)*x(169)-X(50,i)*x(170)-X(51,i)*x(171)-X(52,i)*x(172)-X(53,i)*x(173)-X(54,i)*x(174)-X(55,i)*x(175)-X(56,i)*x(176)-X(57,i)*x(177)-X(58,i)*x(178)-X(59,i)*x(179)-X(60,i)*x(180)-X(61,i)*x(181)-X(62,i)*x(182)-X(63,i)*x(183)-X(64,i)*x(184)-X(65,i)*x(185)-X(66,i)*x(186)-X(67,i)*x(187)-X(68,i)*x(188)-X(69,i)*x(189)-X(70,i)*x(190)-X(71,i)*x(191)-X(72,i)*x(192)-X(73,i)*x(193)-X(74,i)*x(194)-X(75,i)*x(195)-X(76,i)*x(196)-X(77,i)*x(197)-X(78,i)*x(198)-X(79,i)*x(199)-X(80,i)*x(200)-X(81,i)*x(201)-X(82,i)*x(202)-X(83,i)*x(203)-X(84,i)*x(204)-X(85,i)*x(205)-X(86,i)*x(206)-X(87,i)*x(207)-X(88,i)*x(208)-X(89,i)*x(209)-X(90,i)*x(210)-X(91,i)*x(211)-X(92,i)*x(212)-X(93,i)*x(213)-X(94,i)*x(214)-X(95,i)*x(215)-X(96,i)*x(216)-X(97,i)*x(217)-X(98,i)*x(218)-X(99,i)*x(219)-X(100,i)*x(220)-X(101,i)*x(221)-X(102,i)*x(222)-X(103,i)*x(223)-X(104,i)*x(224)-X(105,i)*x(225)-X(106,i)*x(226)-X(107,i)*x(227)-X(108,i)*x(228)-X(109,i)*x(229)-X(110,i)*x(230)-X(111,i)*x(231)-X(112,i)*x(232)-X(113,i)*x(233)-X(114,i)*x(234)-X(115,i)*x(235)-X(116,i)*x(236)-X(117,i)*x(237)-X(118,i)*x(238)-X(119,i)*x(239)-X(120,i)*x(240))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(241)-X(2,i)*x(242)-X(3,i)*x(243)-X(4,i)*x(244)-X(5,i)*x(245)-X(6,i)*x(246)-X(7,i)*x(247)-X(8,i)*x(248)-X(9,i)*x(249)-X(10,i)*x(250)-X(11,i)*x(251)-X(12,i)*x(252)-X(13,i)*x(253)-X(14,i)*x(254)-X(15,i)*x(255)-X(16,i)*x(256)-X(17,i)*x(257)-X(18,i)*x(258)-X(19,i)*x(259)-X(20,i)*x(260)-X(21,i)*x(261)-X(22,i)*x(262)-X(23,i)*x(263)-X(24,i)*x(264)-X(25,i)*x(265)-X(26,i)*x(266)-X(27,i)*x(267)-X(28,i)*x(268)-X(29,i)*x(269)-X(30,i)*x(270)-X(31,i)*x(271)-X(32,i)*x(272)-X(33,i)*x(273)-X(34,i)*x(274)-X(35,i)*x(275)-X(36,i)*x(276)-X(37,i)*x(277)-X(38,i)*x(278)-X(39,i)*x(279)-X(40,i)*x(280)-X(41,i)*x(281)-X(42,i)*x(282)-X(43,i)*x(283)-X(44,i)*x(284)-X(45,i)*x(285)-X(46,i)*x(286)-X(47,i)*x(287)-X(48,i)*x(288)-X(49,i)*x(289)-X(50,i)*x(290)-X(51,i)*x(291)-X(52,i)*x(292)-X(53,i)*x(293)-X(54,i)*x(294)-X(55,i)*x(295)-X(56,i)*x(296)-X(57,i)*x(297)-X(58,i)*x(298)-X(59,i)*x(299)-X(60,i)*x(300)-X(61,i)*x(301)-X(62,i)*x(302)-X(63,i)*x(303)-X(64,i)*x(304)-X(65,i)*x(305)-X(66,i)*x(306)-X(67,i)*x(307)-X(68,i)*x(308)-X(69,i)*x(309)-X(70,i)*x(310)-X(71,i)*x(311)-X(72,i)*x(312)-X(73,i)*x(313)-X(74,i)*x(314)-X(75,i)*x(315)-X(76,i)*x(316)-X(77,i)*x(317)-X(78,i)*x(318)-X(79,i)*x(319)-X(80,i)*x(320)-X(81,i)*x(321)-X(82,i)*x(322)-X(83,i)*x(323)-X(84,i)*x(324)-X(85,i)*x(325)-X(86,i)*x(326)-X(87,i)*x(327)-X(88,i)*x(328)-X(89,i)*x(329)-X(90,i)*x(330)-X(91,i)*x(331)-X(92,i)*x(332)-X(93,i)*x(333)-X(94,i)*x(334)-X(95,i)*x(335)-X(96,i)*x(336)-X(97,i)*x(337)-X(98,i)*x(338)-X(99,i)*x(339)-X(100,i)*x(340)-X(101,i)*x(341)-X(102,i)*x(342)-X(103,i)*x(343)-X(104,i)*x(344)-X(105,i)*x(345)-X(106,i)*x(346)-X(107,i)*x(347)-X(108,i)*x(348)-X(109,i)*x(349)-X(110,i)*x(350)-X(111,i)*x(351)-X(112,i)*x(352)-X(113,i)*x(353)-X(114,i)*x(354)-X(115,i)*x(355)-X(116,i)*x(356)-X(117,i)*x(357)-X(118,i)*x(358)-X(119,i)*x(359)-X(120,i)*x(360))*norm(Nvalues(i)))^2;
end

end

end

