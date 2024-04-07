function [AIntra, fvalIntra] = GetATernaryForth(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';u = double(u);
X = ones(35, size(ROIpositions,1));
X(1,:) = (ROIpositions(:,1).^4);
X(2,:) = (ROIpositions(:,2).^4);
X(3,:) = (ROIpositions(:,3).^4);
X(4,:) = (ROIpositions(:,1).^3).*ROIpositions(:,2);
X(5,:) = (ROIpositions(:,1).^3).*ROIpositions(:,3);
X(6,:) = (ROIpositions(:,2).^3).*ROIpositions(:,1);
X(7,:) = (ROIpositions(:,2).^3).*ROIpositions(:,3);
X(8,:) = (ROIpositions(:,3).^3).*ROIpositions(:,1);
X(9,:) = (ROIpositions(:,3).^3).*ROIpositions(:,2);
X(10,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,2).^2);
X(11,:) = (ROIpositions(:,1).^2).*(ROIpositions(:,3).^2);
X(12,:) = (ROIpositions(:,2).^2).*(ROIpositions(:,3).^2);
X(13,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2).*ROIpositions(:,3);
X(14,:) = ROIpositions(:,1).*(ROIpositions(:,2).^2).*ROIpositions(:,3);
X(15,:) = ROIpositions(:,1).*ROIpositions(:,2).*(ROIpositions(:,3).^2);
X(16,:) = (ROIpositions(:,1).^3);
X(17,:) = (ROIpositions(:,2).^3);
X(18,:) = (ROIpositions(:,3).^3);
X(19,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2);
X(20,:) = (ROIpositions(:,1).^2).*ROIpositions(:,3);
X(21,:) = (ROIpositions(:,2).^2).*ROIpositions(:,1);
X(22,:) = (ROIpositions(:,2).^2).*ROIpositions(:,3);
X(23,:) = (ROIpositions(:,3).^2).*ROIpositions(:,1);
X(24,:) = (ROIpositions(:,3).^2).*ROIpositions(:,2);
X(25,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);
X(26,:) = (ROIpositions(:,1).^2);
X(27,:) = (ROIpositions(:,2).^2);
X(28,:) = (ROIpositions(:,3).^2);
X(29,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(30,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(31,:) = ROIpositions(:,2).*ROIpositions(:,3);
X(32,:) = ROIpositions(:,1);
X(33,:) = ROIpositions(:,2);
X(34,:) = ROIpositions(:,3);

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];
AInitial = double(AInitial);
% Subject to ...
Aeq = zeros(20,105);
Aeq(1,1) = 4; Aeq(1,39) = 1; Aeq(1,75) = 1; 
Aeq(2,6) = 1; Aeq(2,37) = 4; Aeq(2,77) = 1; 
Aeq(3,8) = 1; Aeq(3,44) = 1; Aeq(3,73) = 4; 

Aeq(4,4) = 3; Aeq(4,45) = 2; Aeq(4,83) = 1; 
Aeq(5,5) = 3; Aeq(5,48) = 1; Aeq(5,81) = 2; 
Aeq(6,10) = 2; Aeq(6,43) = 3; Aeq(6,84) = 1; 
Aeq(7,14) = 1; Aeq(7,42) = 3; Aeq(7,82) = 2;
Aeq(8,11) = 2; Aeq(8,50) = 1; Aeq(8,75) = 3; 
Aeq(9,15) = 1; Aeq(9,47) = 2; Aeq(9,79) = 3; 
Aeq(10,13) = 2; Aeq(10,49) = 2; Aeq(10,85) = 2; 

Aeq(11,16) = 3; Aeq(11,54) = 1; Aeq(11,90) = 1;
Aeq(12,21) = 1; Aeq(12,52) = 3; Aeq(12,92) = 1; 
Aeq(13,23) = 1; Aeq(13,59) = 1; Aeq(13,88) = 3; 

Aeq(14,19) = 2; Aeq(14,56) = 2; Aeq(14,95) = 1; 
Aeq(15,20) = 2; Aeq(15,60) = 1; Aeq(15,93) = 2; 
Aeq(16,25) = 2; Aeq(16,57) = 1; Aeq(16,94) = 2; 

Aeq(17,26) = 2; Aeq(17,64) = 1; Aeq(17,100) = 1; 
Aeq(18,29) = 1; Aeq(18,62) = 2; Aeq(18,101) = 1; 
Aeq(19,30) = 1; Aeq(19,66) = 1; Aeq(19,98) = 2; 
Aeq(20,32) = 1; Aeq(20,68) = 1; Aeq(20,104) = 1; 
Beq = zeros(20,1);

% Optimation
options = optimoptions('fmincon','Display','off');
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
AIntra = reshape(AIntra, 35, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20)-X(21,i)*x(21)-X(22,i)*x(22)-X(23,i)*x(23)-X(24,i)*x(24)-X(25,i)*x(25)-X(26,i)*x(26)-X(27,i)*x(27)-X(28,i)*x(28)-X(29,i)*x(29)-X(30,i)*x(30)-X(31,i)*x(31)-X(32,i)*x(32)-X(33,i)*x(33)-X(34,i)*x(34)-X(35,i)*x(35))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(36)-X(2,i)*x(37)-X(3,i)*x(38)-X(4,i)*x(39)-X(5,i)*x(40)-X(6,i)*x(41)-X(7,i)*x(42)-X(8,i)*x(43)-X(9,i)*x(44)-X(10,i)*x(45)-X(11,i)*x(46)-X(12,i)*x(47)-X(13,i)*x(48)-X(14,i)*x(49)-X(15,i)*x(50)-X(16,i)*x(51)-X(17,i)*x(52)-X(18,i)*x(53)-X(19,i)*x(54)-X(20,i)*x(55)-X(21,i)*x(56)-X(22,i)*x(57)-X(23,i)*x(58)-X(24,i)*x(59)-X(25,i)*x(60)-X(26,i)*x(61)-X(27,i)*x(62)-X(28,i)*x(63)-X(29,i)*x(64)-X(30,i)*x(65)-X(31,i)*x(66)-X(32,i)*x(67)-X(33,i)*x(68)-X(34,i)*x(69)-X(35,i)*x(70))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(71)-X(2,i)*x(72)-X(3,i)*x(73)-X(4,i)*x(74)-X(5,i)*x(75)-X(6,i)*x(76)-X(7,i)*x(77)-X(8,i)*x(78)-X(9,i)*x(79)-X(10,i)*x(80)-X(11,i)*x(81)-X(12,i)*x(82)-X(13,i)*x(83)-X(14,i)*x(84)-X(15,i)*x(85)-X(16,i)*x(86)-X(17,i)*x(87)-X(18,i)*x(88)-X(19,i)*x(89)-X(20,i)*x(90)-X(21,i)*x(91)-X(22,i)*x(92)-X(23,i)*x(93)-X(24,i)*x(94)-X(25,i)*x(95)-X(26,i)*x(96)-X(27,i)*x(97)-X(28,i)*x(98)-X(29,i)*x(99)-X(30,i)*x(100)-X(31,i)*x(101)-X(32,i)*x(102)-X(33,i)*x(103)-X(34,i)*x(104)-X(35,i)*x(105))*norm(Nvalues(i)))^2;
end
    
end

end

