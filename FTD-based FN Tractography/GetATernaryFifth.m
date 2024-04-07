function [AIntra, fvalIntra] = GetATernaryFifth(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';
X = ones(56, size(ROIpositions,1));

X(1,:) = (ROIpositions(:,1).^5);
X(2,:) = (ROIpositions(:,2).^5);
X(3,:) = (ROIpositions(:,3).^5);

X(4,:) = (ROIpositions(:,1).^4).*ROIpositions(:,2);
X(5,:) = (ROIpositions(:,1).^4).*ROIpositions(:,3);
X(6,:) = (ROIpositions(:,2).^4).*ROIpositions(:,1);
X(7,:) = (ROIpositions(:,2).^4).*ROIpositions(:,3);
X(8,:) = (ROIpositions(:,3).^4).*ROIpositions(:,1);
X(9,:) = (ROIpositions(:,3).^4).*ROIpositions(:,2);

X(10,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,2).^2);
X(11,:) = (ROIpositions(:,1).^3).*(ROIpositions(:,3).^2);
X(12,:) = (ROIpositions(:,2).^3).*(ROIpositions(:,1).^2);
X(13,:) = (ROIpositions(:,2).^3).*(ROIpositions(:,3).^2);
X(14,:) = (ROIpositions(:,3).^3).*(ROIpositions(:,1).^2);
X(15,:) = (ROIpositions(:,3).^3).*(ROIpositions(:,2).^2);

X(16,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2).*ROIpositions(:,3);
X(17,:) = ROIpositions(:,1).*(ROIpositions(:,2).^3).*ROIpositions(:,3);
X(18,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^3);

X(19,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(20,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(21,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);

X(22,:) = (ROIpositions(:,1).^4);
X(23,:) = (ROIpositions(:,2).^4);
X(24,:) = (ROIpositions(:,3).^4);

X(25,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2);
X(26,:) = (ROIpositions(:,1).^3).*ROIpositions(:,3);
X(27,:) = (ROIpositions(:,2).^3).*ROIpositions(:,1);
X(28,:) = (ROIpositions(:,2).^3).*ROIpositions(:,3);
X(29,:) = (ROIpositions(:,3).^3).*ROIpositions(:,1);
X(30,:) = (ROIpositions(:,3).^3).*ROIpositions(:,2);

X(31,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2);
X(32,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,3).^2);
X(33,:) = (ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);

X(34,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*ROIpositions(:,3);
X(35,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(36,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^2);

X(37,:) = (ROIpositions(:,1).^3);
X(38,:) = (ROIpositions(:,2).^3);
X(39,:) = (ROIpositions(:,3).^3);

X(40,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2);
X(41,:) = (ROIpositions(:,1).^2).*ROIpositions(:,3);
X(42,:) = (ROIpositions(:,2).^2).*ROIpositions(:,1);
X(43,:) = (ROIpositions(:,2).^2).*ROIpositions(:,3);
X(44,:) = (ROIpositions(:,3).^2).*ROIpositions(:,1);
X(45,:) = (ROIpositions(:,3).^2).*ROIpositions(:,2);

X(46,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);

X(47,:) = (ROIpositions(:,1).^2);
X(48,:) = (ROIpositions(:,2).^2);
X(49,:) = (ROIpositions(:,3).^2);

X(50,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(51,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(52,:) = ROIpositions(:,2).*ROIpositions(:,3);

X(53,:) = ROIpositions(:,1);
X(54,:) = ROIpositions(:,2);
X(55,:) = ROIpositions(:,3);

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];

% Subject to ...
Aeq = zeros(35,168);

Aeq(1,1) = 5; Aeq(1,60) = 1; Aeq(1,117) = 1; 
Aeq(2,6) = 1; Aeq(2,58) = 5; Aeq(2,119) = 1; 
Aeq(3,8) = 1; Aeq(3,65) = 1; Aeq(3,115) = 4; 

Aeq(4,4) = 4; Aeq(4,66) = 2; Aeq(4,128) = 1; 
Aeq(5,5) = 4; Aeq(5,72) = 1; Aeq(5,123) = 2; 
Aeq(6,12) = 2; Aeq(6,58) = 4; Aeq(6,129) = 1; 
Aeq(7,17) = 1; Aeq(7,63) = 4; Aeq(7,125) = 2;
Aeq(8,14) = 2; Aeq(8,74) = 1; Aeq(8,120) = 4; 
Aeq(9,18) = 1; Aeq(9,71) = 2; Aeq(9,121) = 4; 

Aeq(10,10) = 3; Aeq(10,58) = 3; Aeq(10,131) = 1; 
Aeq(11,11) = 3; Aeq(11,76) = 1; Aeq(11,126) = 3; 
Aeq(12,21) = 1; Aeq(12,69) = 3; Aeq(12,127) = 3; 

Aeq(13,16) = 3; Aeq(13,75) = 2; Aeq(13,132) = 2;
Aeq(14,19) = 2; Aeq(14,73) = 3; Aeq(14,133) = 2; 
Aeq(15,20) = 2; Aeq(15,77) = 2; Aeq(15,130) = 3; 

Aeq(16,22) = 4; Aeq(16,81) = 1; Aeq(16,138) = 1; 
Aeq(17,27) = 1; Aeq(17,79) = 4; Aeq(17,140) = 1; 
Aeq(18,29) = 1; Aeq(18,86) = 1; Aeq(18,136) = 4; 

Aeq(19,25) = 3; Aeq(19,87) = 2; Aeq(19,146) = 1; 
Aeq(20,26) = 3; Aeq(20,90) = 1; Aeq(20,144) = 2; 
Aeq(21,31) = 2; Aeq(21,85) = 3; Aeq(21,147) = 1; 
Aeq(22,35) = 1; Aeq(22,84) = 3; Aeq(22,145) = 2;
Aeq(23,32) = 2; Aeq(23,92) = 1; Aeq(23,138) = 3; 
Aeq(24,36) = 1; Aeq(24,89) = 2; Aeq(24,142) = 3; 
Aeq(25,34) = 2; Aeq(25,91) = 2; Aeq(25,148) = 2; 

Aeq(26,37) = 3; Aeq(26,96) = 1; Aeq(26,153) = 1;
Aeq(27,42) = 1; Aeq(27,94) = 3; Aeq(27,155) = 1; 
Aeq(28,44) = 1; Aeq(28,101) = 1; Aeq(28,151) = 3; 

Aeq(29,40) = 2; Aeq(29,98) = 2; Aeq(29,158) = 1; 
Aeq(30,41) = 2; Aeq(30,102) = 1; Aeq(30,156) = 2; 
Aeq(31,46) = 2; Aeq(31,99) = 1; Aeq(31,157) = 2; 

Aeq(32,47) = 2; Aeq(32,106) = 1; Aeq(32,163) = 1; 
Aeq(33,50) = 1; Aeq(33,104) = 2; Aeq(33,164) = 1; 
Aeq(34,51) = 1; Aeq(34,108) = 1; Aeq(34,161) = 2; 
Aeq(35,53) = 1; Aeq(35,110) = 1; Aeq(35,167) = 1; 

Beq = zeros(35,1);

% Optimation
options = optimoptions('fmincon','Display','off');
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
for i=1:size(X,2)
    a = -X(1,i)*AIntra(1)-X(2,i)*AIntra(2)-X(3,i)*AIntra(3)-X(4,i)*AIntra(4)-X(5,i)*AIntra(5)-X(6,i)*AIntra(6)-X(7,i)*AIntra(7)-X(8,i)*AIntra(8)-X(9,i)*AIntra(9)-X(10,i)*AIntra(10)-X(11,i)*AIntra(11)-X(12,i)*AIntra(12)-X(13,i)*AIntra(13)-X(14,i)*AIntra(14)-X(15,i)*AIntra(15)-X(16,i)*AIntra(16)-X(17,i)*AIntra(17)-X(18,i)*AIntra(18)-X(19,i)*AIntra(19)-X(20,i)*AIntra(20)-X(21,i)*AIntra(21)-X(22,i)*AIntra(22)-X(23,i)*AIntra(23)-X(24,i)*AIntra(24)-X(25,i)*AIntra(25)-X(26,i)*AIntra(26)-X(27,i)*AIntra(27)-X(28,i)*AIntra(28)-X(29,i)*AIntra(29)-X(30,i)*AIntra(30)-X(31,i)*AIntra(31)-X(32,i)*AIntra(32)-X(33,i)*AIntra(33)-X(34,i)*AIntra(34)-X(35,i)*AIntra(35)-X(36,i)*AIntra(36)-X(37,i)*AIntra(37)-X(38,i)*AIntra(38)-X(39,i)*AIntra(39)-X(40,i)*AIntra(40)-X(41,i)*AIntra(41)-X(42,i)*AIntra(42)-X(43,i)*AIntra(43)-X(44,i)*AIntra(44)-X(45,i)*AIntra(45)-X(46,i)*AIntra(46)-X(47,i)*AIntra(47)-X(48,i)*AIntra(48)-X(49,i)*AIntra(49)-X(50,i)*AIntra(50)-X(51,i)*AIntra(51)-X(52,i)*AIntra(52)-X(53,i)*AIntra(53)-X(54,i)*AIntra(54)-X(55,i)*AIntra(55)-X(56,i)*AIntra(56);
end
AIntra = reshape(AIntra, 56, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20)-X(21,i)*x(21)-X(22,i)*x(22)-X(23,i)*x(23)-X(24,i)*x(24)-X(25,i)*x(25)-X(26,i)*x(26)-X(27,i)*x(27)-X(28,i)*x(28)-X(29,i)*x(29)-X(30,i)*x(30)-X(31,i)*x(31)-X(32,i)*x(32)-X(33,i)*x(33)-X(34,i)*x(34)-X(35,i)*x(35)-X(36,i)*x(36)-X(37,i)*x(37)-X(38,i)*x(38)-X(39,i)*x(39)-X(40,i)*x(40)-X(41,i)*x(41)-X(42,i)*x(42)-X(43,i)*x(43)-X(44,i)*x(44)-X(45,i)*x(45)-X(46,i)*x(46)-X(47,i)*x(47)-X(48,i)*x(48)-X(49,i)*x(49)-X(50,i)*x(50)-X(51,i)*x(51)-X(52,i)*x(52)-X(53,i)*x(53)-X(54,i)*x(54)-X(55,i)*x(55)-X(56,i)*x(56))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(57)-X(2,i)*x(58)-X(3,i)*x(59)-X(4,i)*x(60)-X(5,i)*x(61)-X(6,i)*x(62)-X(7,i)*x(63)-X(8,i)*x(64)-X(9,i)*x(65)-X(10,i)*x(66)-X(11,i)*x(67)-X(12,i)*x(68)-X(13,i)*x(69)-X(14,i)*x(70)-X(15,i)*x(71)-X(16,i)*x(72)-X(17,i)*x(73)-X(18,i)*x(74)-X(19,i)*x(75)-X(20,i)*x(76)-X(21,i)*x(77)-X(22,i)*x(78)-X(23,i)*x(79)-X(24,i)*x(80)-X(25,i)*x(81)-X(26,i)*x(82)-X(27,i)*x(83)-X(28,i)*x(84)-X(29,i)*x(85)-X(30,i)*x(86)-X(31,i)*x(87)-X(32,i)*x(88)-X(33,i)*x(89)-X(34,i)*x(90)-X(35,i)*x(91)-X(36,i)*x(92)-X(37,i)*x(93)-X(38,i)*x(94)-X(39,i)*x(95)-X(40,i)*x(96)-X(41,i)*x(97)-X(42,i)*x(98)-X(43,i)*x(99)-X(44,i)*x(100)-X(45,i)*x(101)-X(46,i)*x(102)-X(47,i)*x(103)-X(48,i)*x(104)-X(49,i)*x(105)-X(50,i)*x(106)-X(51,i)*x(107)-X(52,i)*x(108)-X(53,i)*x(109)-X(54,i)*x(110)-X(55,i)*x(111)-X(56,i)*x(112))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(113)-X(2,i)*x(114)-X(3,i)*x(115)-X(4,i)*x(116)-X(5,i)*x(117)-X(6,i)*x(118)-X(7,i)*x(119)-X(8,i)*x(120)-X(9,i)*x(121)-X(10,i)*x(122)-X(11,i)*x(123)-X(12,i)*x(124)-X(13,i)*x(125)-X(14,i)*x(126)-X(15,i)*x(127)-X(16,i)*x(128)-X(17,i)*x(129)-X(18,i)*x(130)-X(19,i)*x(131)-X(20,i)*x(132)-X(21,i)*x(133)-X(22,i)*x(134)-X(23,i)*x(135)-X(24,i)*x(136)-X(25,i)*x(137)-X(26,i)*x(138)-X(27,i)*x(139)-X(28,i)*x(140)-X(29,i)*x(141)-X(30,i)*x(142)-X(31,i)*x(143)-X(32,i)*x(144)-X(33,i)*x(145)-X(34,i)*x(146)-X(35,i)*x(147)-X(36,i)*x(148)-X(37,i)*x(149)-X(38,i)*x(150)-X(39,i)*x(151)-X(40,i)*x(152)-X(41,i)*x(153)-X(42,i)*x(154)-X(43,i)*x(155)-X(44,i)*x(156)-X(45,i)*x(157)-X(46,i)*x(158)-X(47,i)*x(159)-X(48,i)*x(160)-X(49,i)*x(161)-X(50,i)*x(162)-X(51,i)*x(163)-X(52,i)*x(164)-X(53,i)*x(165)-X(54,i)*x(166)-X(55,i)*x(167)-X(56,i)*x(168))*norm(Nvalues(i)))^2;
end
    
end

end

