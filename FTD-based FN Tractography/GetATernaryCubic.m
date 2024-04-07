function [AIntra, fvalIntra] = GetATernaryCubic(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';
X = ones(20, size(ROIpositions,1));
X(1,:) = (ROIpositions(:,1).^3);
X(2,:) = (ROIpositions(:,2).^3);
X(3,:) = (ROIpositions(:,3).^3);
X(4,:) = (ROIpositions(:,1).^2).*ROIpositions(:,2);
X(5,:) = (ROIpositions(:,1).^2).*ROIpositions(:,3);
X(6,:) = (ROIpositions(:,2).^2).*ROIpositions(:,1);
X(7,:) = (ROIpositions(:,2).^2).*ROIpositions(:,3);
X(8,:) = (ROIpositions(:,3).^2).*ROIpositions(:,1);
X(9,:) = (ROIpositions(:,3).^2).*ROIpositions(:,2);
X(10,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);
X(11,:) = (ROIpositions(:,1).^2);
X(12,:) = (ROIpositions(:,2).^2);
X(13,:) = (ROIpositions(:,3).^2);
X(14,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(15,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(16,:) = ROIpositions(:,2).*ROIpositions(:,3);
X(17,:) = ROIpositions(:,1);
X(18,:) = ROIpositions(:,2);
X(19,:) = ROIpositions(:,3);

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];

% Subject to ...
Aeq = zeros(10,60);
Aeq(1,1) = 3; Aeq(1,24) = 1; Aeq(1,45) = 1; 
Aeq(2,6) = 1; Aeq(2,22) = 3; Aeq(2,47) = 1; 
Aeq(3,8) = 1; Aeq(3,29) = 1; Aeq(3,43) = 3; 
Aeq(4,4) = 2; Aeq(4,26) = 2; Aeq(4,50) = 1; 
Aeq(5,5) = 2; Aeq(5,30) = 1; Aeq(5,48) = 2; 
Aeq(6,10) = 1; Aeq(6,27) = 2; Aeq(6,49) = 2; 
Aeq(7,11) = 2; Aeq(7,34) = 1; Aeq(7,55) = 1; 
Aeq(8,14) = 1; Aeq(8,32) = 2; Aeq(8,56) = 1; 
Aeq(9,15) = 1; Aeq(9,36) = 1; Aeq(9,53) = 2; 
Aeq(10,17) = 1; Aeq(10,38) = 1; Aeq(10,59) = 1; 
Beq = zeros(10,1);

% Optimation
options = optimoptions('fmincon','Display','off');
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
AIntra = reshape(AIntra, 20, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(21)-X(2,i)*x(22)-X(3,i)*x(23)-X(4,i)*x(24)-X(5,i)*x(25)-X(6,i)*x(26)-X(7,i)*x(27)-X(8,i)*x(28)-X(9,i)*x(29)-X(10,i)*x(30)-X(11,i)*x(31)-X(12,i)*x(32)-X(13,i)*x(33)-X(14,i)*x(34)-X(15,i)*x(35)-X(16,i)*x(36)-X(17,i)*x(37)-X(18,i)*x(38)-X(19,i)*x(39)-X(20,i)*x(40))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(41)-X(2,i)*x(42)-X(3,i)*x(43)-X(4,i)*x(44)-X(5,i)*x(45)-X(6,i)*x(46)-X(7,i)*x(47)-X(8,i)*x(48)-X(9,i)*x(49)-X(10,i)*x(50)-X(11,i)*x(51)-X(12,i)*x(52)-X(13,i)*x(53)-X(14,i)*x(54)-X(15,i)*x(55)-X(16,i)*x(56)-X(17,i)*x(57)-X(18,i)*x(58)-X(19,i)*x(59)-X(20,i)*x(60))*norm(Nvalues(i)))^2;
end
    
end

end

