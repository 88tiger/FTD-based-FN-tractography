function [AIntra, fvalIntra] = GetATernaryEighth(ROIpositions, DirsROI, WeightedDirsROI)

DirsROI = double(DirsROI);
WeightedDirsROI = double(WeightedDirsROI);
% Weight
Nvalues = sqrt(sum(WeightedDirsROI.^2,2));

% Define Items x and Scalar-Valued u
u = DirsROI';
X = ones(165, size(ROIpositions,1));

X(1,:) = ROIpositions(:,1);
X(2,:) = ROIpositions(:,1).^2;
X(3,:) = ROIpositions(:,1).^3;
X(4,:) = ROIpositions(:,1).^4;
X(5,:) = ROIpositions(:,1).^5;
X(6,:) = ROIpositions(:,1).^6;
X(7,:) = ROIpositions(:,1).^7;
X(8,:) = ROIpositions(:,1).^8;
X(9,:) = ROIpositions(:,2);
X(10,:) = ROIpositions(:,1).*ROIpositions(:,2);
X(11,:) = ROIpositions(:,1).^2.*ROIpositions(:,2);
X(12,:) = ROIpositions(:,1).^3.*ROIpositions(:,2);
X(13,:) = ROIpositions(:,1).^4.*ROIpositions(:,2);
X(14,:) = ROIpositions(:,1).^5.*ROIpositions(:,2);
X(15,:) = ROIpositions(:,1).^6.*ROIpositions(:,2);
X(16,:) = ROIpositions(:,1).^7.*ROIpositions(:,2);
X(17,:) = ROIpositions(:,2).^2;
X(18,:) = ROIpositions(:,1).*ROIpositions(:,2).^2;
X(19,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2;
X(20,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2;
X(21,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^2;
X(22,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).^2;
X(23,:) = ROIpositions(:,1).^6.*ROIpositions(:,2).^2;
X(24,:) = ROIpositions(:,2).^3;
X(25,:) = ROIpositions(:,1).*ROIpositions(:,2).^3;
X(26,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3;
X(27,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^3;
X(28,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^3;
X(29,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).^3;
X(30,:) = ROIpositions(:,2).^4;
X(31,:) = ROIpositions(:,1).*ROIpositions(:,2).^4;
X(32,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^4;
X(33,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^4;
X(34,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^4;
X(35,:) = ROIpositions(:,2).^5;
X(36,:) = ROIpositions(:,1).*ROIpositions(:,2).^5;
X(37,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^5;
X(38,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^5;
X(39,:) = ROIpositions(:,2).^6;
X(40,:) = ROIpositions(:,1).*ROIpositions(:,2).^6;
X(41,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^6;
X(42,:) = ROIpositions(:,2).^7;
X(43,:) = ROIpositions(:,1).*ROIpositions(:,2).^7;
X(44,:) = ROIpositions(:,2).^8;
X(45,:) = ROIpositions(:,3);
X(46,:) = ROIpositions(:,1).*ROIpositions(:,3);
X(47,:) = ROIpositions(:,1).^2.*ROIpositions(:,3);
X(48,:) = ROIpositions(:,1).^3.*ROIpositions(:,3);
X(49,:) = ROIpositions(:,1).^4.*ROIpositions(:,3);
X(50,:) = ROIpositions(:,1).^5.*ROIpositions(:,3);
X(51,:) = ROIpositions(:,1).^6.*ROIpositions(:,3);
X(52,:) = ROIpositions(:,1).^7.*ROIpositions(:,3);
X(53,:) = ROIpositions(:,2).*ROIpositions(:,3);
X(54,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3);
X(55,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3);
X(56,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3);
X(57,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).*ROIpositions(:,3);
X(58,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).*ROIpositions(:,3);
X(59,:) = ROIpositions(:,1).^6.*ROIpositions(:,2).*ROIpositions(:,3);
X(60,:) = ROIpositions(:,2).^2.*ROIpositions(:,3);
X(61,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(62,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(63,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(64,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(65,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).^2.*ROIpositions(:,3);
X(66,:) = ROIpositions(:,2).^3.*ROIpositions(:,3);
X(67,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(68,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(69,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(70,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^3.*ROIpositions(:,3);
X(71,:) = ROIpositions(:,2).^4.*ROIpositions(:,3);
X(72,:) = ROIpositions(:,1).*ROIpositions(:,2).^4.*ROIpositions(:,3);
X(73,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^4.*ROIpositions(:,3);
X(74,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^4.*ROIpositions(:,3);
X(75,:) = ROIpositions(:,2).^5.*ROIpositions(:,3);
X(76,:) = ROIpositions(:,1).*ROIpositions(:,2).^5.*ROIpositions(:,3);
X(77,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^5.*ROIpositions(:,3);
X(78,:) = ROIpositions(:,2).^6.*ROIpositions(:,3);
X(79,:) = ROIpositions(:,1).*ROIpositions(:,2).^6.*ROIpositions(:,3);
X(80,:) = ROIpositions(:,2).^7.*ROIpositions(:,3);
X(81,:) = ROIpositions(:,3).^2;
X(82,:) = ROIpositions(:,1).*ROIpositions(:,3).^2;
X(83,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^2;
X(84,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^2;
X(85,:) = ROIpositions(:,1).^4.*ROIpositions(:,3).^2;
X(86,:) = ROIpositions(:,1).^5.*ROIpositions(:,3).^2;
X(87,:) = ROIpositions(:,1).^6.*ROIpositions(:,3).^2;
X(88,:) = ROIpositions(:,2).*ROIpositions(:,3).^2;
X(89,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(90,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(91,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(92,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(93,:) = ROIpositions(:,1).^5.*ROIpositions(:,2).*ROIpositions(:,3).^2;
X(94,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(95,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(96,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(97,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(98,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).^2.*ROIpositions(:,3).^2;
X(99,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(100,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(101,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(102,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^3.*ROIpositions(:,3).^2;
X(103,:) = ROIpositions(:,2).^4.*ROIpositions(:,3).^2;
X(104,:) = ROIpositions(:,1).*ROIpositions(:,2).^4.*ROIpositions(:,3).^2;
X(105,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^4.*ROIpositions(:,3).^2;
X(106,:) = ROIpositions(:,2).^5.*ROIpositions(:,3).^2;
X(107,:) = ROIpositions(:,1).*ROIpositions(:,2).^5.*ROIpositions(:,3).^2;
X(108,:) = ROIpositions(:,2).^6.*ROIpositions(:,3).^2;
X(109,:) = ROIpositions(:,3).^3;
X(110,:) = ROIpositions(:,1).*ROIpositions(:,3).^3;
X(111,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^3;
X(112,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^3;
X(113,:) = ROIpositions(:,1).^4.*ROIpositions(:,3).^3;
X(114,:) = ROIpositions(:,1).^5.*ROIpositions(:,3).^3;
X(115,:) = ROIpositions(:,2).*ROIpositions(:,3).^3;
X(116,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(117,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(118,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(119,:) = ROIpositions(:,1).^4.*ROIpositions(:,2).*ROIpositions(:,3).^3;
X(120,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(121,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(122,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(123,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).^2.*ROIpositions(:,3).^3;
X(124,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^3;
X(125,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3).^3;
X(126,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^3.*ROIpositions(:,3).^3;
X(127,:) = ROIpositions(:,2).^4.*ROIpositions(:,3).^3;
X(128,:) = ROIpositions(:,1).*ROIpositions(:,2).^4.*ROIpositions(:,3).^3;
X(129,:) = ROIpositions(:,2).^5.*ROIpositions(:,3).^3;
X(130,:) = ROIpositions(:,3).^4;
X(131,:) = ROIpositions(:,1).*ROIpositions(:,3).^4;
X(132,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^4;
X(133,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^4;
X(134,:) = ROIpositions(:,1).^4.*ROIpositions(:,3).^4;
X(135,:) = ROIpositions(:,2).*ROIpositions(:,3).^4;
X(136,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^4;
X(137,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^4;
X(138,:) = ROIpositions(:,1).^3.*ROIpositions(:,2).*ROIpositions(:,3).^4;
X(139,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^4;
X(140,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^4;
X(141,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).^2.*ROIpositions(:,3).^4;
X(142,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^4;
X(143,:) = ROIpositions(:,1).*ROIpositions(:,2).^3.*ROIpositions(:,3).^4;
X(144,:) = ROIpositions(:,2).^4.*ROIpositions(:,3).^4;
X(145,:) = ROIpositions(:,3).^5;
X(146,:) = ROIpositions(:,1).*ROIpositions(:,3).^5;
X(147,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^5;
X(148,:) = ROIpositions(:,1).^3.*ROIpositions(:,3).^5;
X(149,:) = ROIpositions(:,2).*ROIpositions(:,3).^5;
X(150,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^5;
X(151,:) = ROIpositions(:,1).^2.*ROIpositions(:,2).*ROIpositions(:,3).^5;
X(152,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^5;
X(153,:) = ROIpositions(:,1).*ROIpositions(:,2).^2.*ROIpositions(:,3).^5;
X(154,:) = ROIpositions(:,2).^3.*ROIpositions(:,3).^5;
X(155,:) = ROIpositions(:,3).^6;
X(156,:) = ROIpositions(:,1).*ROIpositions(:,3).^6;
X(157,:) = ROIpositions(:,1).^2.*ROIpositions(:,3).^6;
X(158,:) = ROIpositions(:,2).*ROIpositions(:,3).^6;
X(159,:) = ROIpositions(:,1).*ROIpositions(:,2).*ROIpositions(:,3).^6;
X(160,:) = ROIpositions(:,2).^2.*ROIpositions(:,3).^6;
X(161,:) = ROIpositions(:,3).^7;
X(162,:) = ROIpositions(:,1).*ROIpositions(:,3).^7;
X(163,:) = ROIpositions(:,2).*ROIpositions(:,3).^7;
X(164,:) = ROIpositions(:,3).^8;

% Initial A
AInversion = u*pinv(X);
AInitial = [AInversion(1,:), AInversion(2,:), AInversion(3,:)];

% Subject to ...
Aeq = zeros(120,495);

Aeq(1,1) = 1; Aeq(1,174) = 1; Aeq(1,375) = 1;
Aeq(2,2) = 2; Aeq(2,175) = 1; Aeq(2,376) = 1;
Aeq(3,3) = 3; Aeq(3,176) = 1; Aeq(3,377) = 1;
Aeq(4,4) = 4; Aeq(4,177) = 1; Aeq(4,378) = 1;
Aeq(5,5) = 5; Aeq(5,178) = 1; Aeq(5,379) = 1;
Aeq(6,6) = 6; Aeq(6,179) = 1; Aeq(6,380) = 1;
Aeq(7,7) = 7; Aeq(7,180) = 1; Aeq(7,381) = 1;
Aeq(8,8) = 8; Aeq(8,181) = 1; Aeq(8,382) = 1;
Aeq(9,10) = 1; Aeq(9,182) = 2; Aeq(9,383) = 1;
Aeq(10,11) = 2; Aeq(10,183) = 2; Aeq(10,384) = 1;
Aeq(11,12) = 3; Aeq(11,184) = 2; Aeq(11,385) = 1;
Aeq(12,13) = 4; Aeq(12,185) = 2; Aeq(12,386) = 1;
Aeq(13,14) = 5; Aeq(13,186) = 2; Aeq(13,387) = 1;
Aeq(14,15) = 6; Aeq(14,187) = 2; Aeq(14,388) = 1;
Aeq(15,16) = 7; Aeq(15,188) = 2; Aeq(15,389) = 1;
Aeq(16,18) = 1; Aeq(16,189) = 3; Aeq(16,390) = 1;
Aeq(17,19) = 2; Aeq(17,190) = 3; Aeq(17,391) = 1;
Aeq(18,20) = 3; Aeq(18,191) = 3; Aeq(18,392) = 1;
Aeq(19,21) = 4; Aeq(19,192) = 3; Aeq(19,393) = 1;
Aeq(20,22) = 5; Aeq(20,193) = 3; Aeq(20,394) = 1;
Aeq(21,23) = 6; Aeq(21,194) = 3; Aeq(21,395) = 1;
Aeq(22,25) = 1; Aeq(22,195) = 4; Aeq(22,396) = 1;
Aeq(23,26) = 2; Aeq(23,196) = 4; Aeq(23,397) = 1;
Aeq(24,27) = 3; Aeq(24,197) = 4; Aeq(24,398) = 1;
Aeq(25,28) = 4; Aeq(25,198) = 4; Aeq(25,399) = 1;
Aeq(26,29) = 5; Aeq(26,199) = 4; Aeq(26,400) = 1;
Aeq(27,31) = 1; Aeq(27,200) = 5; Aeq(27,401) = 1;
Aeq(28,32) = 2; Aeq(28,201) = 5; Aeq(28,402) = 1;
Aeq(29,33) = 3; Aeq(29,202) = 5; Aeq(29,403) = 1;
Aeq(30,34) = 4; Aeq(30,203) = 5; Aeq(30,404) = 1;
Aeq(31,36) = 1; Aeq(31,204) = 6; Aeq(31,405) = 1;
Aeq(32,37) = 2; Aeq(32,205) = 6; Aeq(32,406) = 1;
Aeq(33,38) = 3; Aeq(33,206) = 6; Aeq(33,407) = 1;
Aeq(34,40) = 1; Aeq(34,207) = 7; Aeq(34,408) = 1;
Aeq(35,41) = 2; Aeq(35,208) = 7; Aeq(35,409) = 1;
Aeq(36,43) = 1; Aeq(36,209) = 8; Aeq(36,410) = 1;
Aeq(37,46) = 1; Aeq(37,218) = 1; Aeq(37,411) = 2;
Aeq(38,47) = 2; Aeq(38,219) = 1; Aeq(38,412) = 2;
Aeq(39,48) = 3; Aeq(39,220) = 1; Aeq(39,413) = 2;
Aeq(40,49) = 4; Aeq(40,221) = 1; Aeq(40,414) = 2;
Aeq(41,50) = 5; Aeq(41,222) = 1; Aeq(41,415) = 2;
Aeq(42,51) = 6; Aeq(42,223) = 1; Aeq(42,416) = 2;
Aeq(43,52) = 7; Aeq(43,224) = 1; Aeq(43,417) = 2;
Aeq(44,54) = 1; Aeq(44,225) = 2; Aeq(44,418) = 2;
Aeq(45,55) = 2; Aeq(45,226) = 2; Aeq(45,419) = 2;
Aeq(46,56) = 3; Aeq(46,227) = 2; Aeq(46,420) = 2;
Aeq(47,57) = 4; Aeq(47,228) = 2; Aeq(47,421) = 2;
Aeq(48,58) = 5; Aeq(48,229) = 2; Aeq(48,422) = 2;
Aeq(49,59) = 6; Aeq(49,230) = 2; Aeq(49,423) = 2;
Aeq(50,61) = 1; Aeq(50,231) = 3; Aeq(50,424) = 2;
Aeq(51,62) = 2; Aeq(51,232) = 3; Aeq(51,425) = 2;
Aeq(52,63) = 3; Aeq(52,233) = 3; Aeq(52,426) = 2;
Aeq(53,64) = 4; Aeq(53,234) = 3; Aeq(53,427) = 2;
Aeq(54,65) = 5; Aeq(54,235) = 3; Aeq(54,428) = 2;
Aeq(55,67) = 1; Aeq(55,236) = 4; Aeq(55,429) = 2;
Aeq(56,68) = 2; Aeq(56,237) = 4; Aeq(56,430) = 2;
Aeq(57,69) = 3; Aeq(57,238) = 4; Aeq(57,431) = 2;
Aeq(58,70) = 4; Aeq(58,239) = 4; Aeq(58,432) = 2;
Aeq(59,72) = 1; Aeq(59,240) = 5; Aeq(59,433) = 2;
Aeq(60,73) = 2; Aeq(60,241) = 5; Aeq(60,434) = 2;
Aeq(61,74) = 3; Aeq(61,242) = 5; Aeq(61,435) = 2;
Aeq(62,76) = 1; Aeq(62,243) = 6; Aeq(62,436) = 2;
Aeq(63,77) = 2; Aeq(63,244) = 6; Aeq(63,437) = 2;
Aeq(64,79) = 1; Aeq(64,245) = 7; Aeq(64,438) = 2;
Aeq(65,82) = 1; Aeq(65,253) = 1; Aeq(65,439) = 3;
Aeq(66,83) = 2; Aeq(66,254) = 1; Aeq(66,440) = 3;
Aeq(67,84) = 3; Aeq(67,255) = 1; Aeq(67,441) = 3;
Aeq(68,85) = 4; Aeq(68,256) = 1; Aeq(68,442) = 3;
Aeq(69,86) = 5; Aeq(69,257) = 1; Aeq(69,443) = 3;
Aeq(70,87) = 6; Aeq(70,258) = 1; Aeq(70,444) = 3;
Aeq(71,89) = 1; Aeq(71,259) = 2; Aeq(71,445) = 3;
Aeq(72,90) = 2; Aeq(72,260) = 2; Aeq(72,446) = 3;
Aeq(73,91) = 3; Aeq(73,261) = 2; Aeq(73,447) = 3;
Aeq(74,92) = 4; Aeq(74,262) = 2; Aeq(74,448) = 3;
Aeq(75,93) = 5; Aeq(75,263) = 2; Aeq(75,449) = 3;
Aeq(76,95) = 1; Aeq(76,264) = 3; Aeq(76,450) = 3;
Aeq(77,96) = 2; Aeq(77,265) = 3; Aeq(77,451) = 3;
Aeq(78,97) = 3; Aeq(78,266) = 3; Aeq(78,452) = 3;
Aeq(79,98) = 4; Aeq(79,267) = 3; Aeq(79,453) = 3;
Aeq(80,100) = 1; Aeq(80,268) = 4; Aeq(80,454) = 3;
Aeq(81,101) = 2; Aeq(81,269) = 4; Aeq(81,455) = 3;
Aeq(82,102) = 3; Aeq(82,270) = 4; Aeq(82,456) = 3;
Aeq(83,104) = 1; Aeq(83,271) = 5; Aeq(83,457) = 3;
Aeq(84,105) = 2; Aeq(84,272) = 5; Aeq(84,458) = 3;
Aeq(85,107) = 1; Aeq(85,273) = 6; Aeq(85,459) = 3;
Aeq(86,110) = 1; Aeq(86,280) = 1; Aeq(86,460) = 4;
Aeq(87,111) = 2; Aeq(87,281) = 1; Aeq(87,461) = 4;
Aeq(88,112) = 3; Aeq(88,282) = 1; Aeq(88,462) = 4;
Aeq(89,113) = 4; Aeq(89,283) = 1; Aeq(89,463) = 4;
Aeq(90,114) = 5; Aeq(90,284) = 1; Aeq(90,464) = 4;
Aeq(91,116) = 1; Aeq(91,285) = 2; Aeq(91,465) = 4;
Aeq(92,117) = 2; Aeq(92,286) = 2; Aeq(92,466) = 4;
Aeq(93,118) = 3; Aeq(93,287) = 2; Aeq(93,467) = 4;
Aeq(94,119) = 4; Aeq(94,288) = 2; Aeq(94,468) = 4;
Aeq(95,121) = 1; Aeq(95,289) = 3; Aeq(95,469) = 4;
Aeq(96,122) = 2; Aeq(96,290) = 3; Aeq(96,470) = 4;
Aeq(97,123) = 3; Aeq(97,291) = 3; Aeq(97,471) = 4;
Aeq(98,125) = 1; Aeq(98,292) = 4; Aeq(98,472) = 4;
Aeq(99,126) = 2; Aeq(99,293) = 4; Aeq(99,473) = 4;
Aeq(100,128) = 1; Aeq(100,294) = 5; Aeq(100,474) = 4;
Aeq(101,131) = 1; Aeq(101,300) = 1; Aeq(101,475) = 5;
Aeq(102,132) = 2; Aeq(102,301) = 1; Aeq(102,476) = 5;
Aeq(103,133) = 3; Aeq(103,302) = 1; Aeq(103,477) = 5;
Aeq(104,134) = 4; Aeq(104,303) = 1; Aeq(104,478) = 5;
Aeq(105,136) = 1; Aeq(105,304) = 2; Aeq(105,479) = 5;
Aeq(106,137) = 2; Aeq(106,305) = 2; Aeq(106,480) = 5;
Aeq(107,138) = 3; Aeq(107,306) = 2; Aeq(107,481) = 5;
Aeq(108,140) = 1; Aeq(108,307) = 3; Aeq(108,482) = 5;
Aeq(109,141) = 2; Aeq(109,308) = 3; Aeq(109,483) = 5;
Aeq(110,143) = 1; Aeq(110,309) = 4; Aeq(110,484) = 5;
Aeq(111,146) = 1; Aeq(111,314) = 1; Aeq(111,485) = 6;
Aeq(112,147) = 2; Aeq(112,315) = 1; Aeq(112,486) = 6;
Aeq(113,148) = 3; Aeq(113,316) = 1; Aeq(113,487) = 6;
Aeq(114,150) = 1; Aeq(114,317) = 2; Aeq(114,488) = 6;
Aeq(115,151) = 2; Aeq(115,318) = 2; Aeq(115,489) = 6;
Aeq(116,153) = 1; Aeq(116,319) = 3; Aeq(116,490) = 6;
Aeq(117,156) = 1; Aeq(117,323) = 1; Aeq(117,491) = 7;
Aeq(118,157) = 2; Aeq(118,324) = 1; Aeq(118,492) = 7;
Aeq(119,159) = 1; Aeq(119,325) = 2; Aeq(119,493) = 7;
Aeq(120,162) = 1; Aeq(120,328) = 1; Aeq(120,494) = 8;

Beq = zeros(120,1);

% Optimation
options = optimoptions('fmincon','Display','off');
[AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], Aeq, Beq, [], [], [], options);
% [AIntra, fvalIntra] = fmincon(@fun, AInitial, [], [], [], [], [], [], [], options);
AIntra = reshape(AIntra, 165, 3)';

% Optimation function
function f = fun(x)
    
f = 0;
for i = 1:size(X,2)
    f = f + ((u(1,i)-X(1,i)*x(1)-X(2,i)*x(2)-X(3,i)*x(3)-X(4,i)*x(4)-X(5,i)*x(5)-X(6,i)*x(6)-X(7,i)*x(7)-X(8,i)*x(8)-X(9,i)*x(9)-X(10,i)*x(10)-X(11,i)*x(11)-X(12,i)*x(12)-X(13,i)*x(13)-X(14,i)*x(14)-X(15,i)*x(15)-X(16,i)*x(16)-X(17,i)*x(17)-X(18,i)*x(18)-X(19,i)*x(19)-X(20,i)*x(20)-X(21,i)*x(21)-X(22,i)*x(22)-X(23,i)*x(23)-X(24,i)*x(24)-X(25,i)*x(25)-X(26,i)*x(26)-X(27,i)*x(27)-X(28,i)*x(28)-X(29,i)*x(29)-X(30,i)*x(30)-X(31,i)*x(31)-X(32,i)*x(32)-X(33,i)*x(33)-X(34,i)*x(34)-X(35,i)*x(35)-X(36,i)*x(36)-X(37,i)*x(37)-X(38,i)*x(38)-X(39,i)*x(39)-X(40,i)*x(40)-X(41,i)*x(41)-X(42,i)*x(42)-X(43,i)*x(43)-X(44,i)*x(44)-X(45,i)*x(45)-X(46,i)*x(46)-X(47,i)*x(47)-X(48,i)*x(48)-X(49,i)*x(49)-X(50,i)*x(50)-X(51,i)*x(51)-X(52,i)*x(52)-X(53,i)*x(53)-X(54,i)*x(54)-X(55,i)*x(55)-X(56,i)*x(56)-X(57,i)*x(57)-X(58,i)*x(58)-X(59,i)*x(59)-X(60,i)*x(60)-X(61,i)*x(61)-X(62,i)*x(62)-X(63,i)*x(63)-X(64,i)*x(64)-X(65,i)*x(65)-X(66,i)*x(66)-X(67,i)*x(67)-X(68,i)*x(68)-X(69,i)*x(69)-X(70,i)*x(70)-X(71,i)*x(71)-X(72,i)*x(72)-X(73,i)*x(73)-X(74,i)*x(74)-X(75,i)*x(75)-X(76,i)*x(76)-X(77,i)*x(77)-X(78,i)*x(78)-X(79,i)*x(79)-X(80,i)*x(80)-X(81,i)*x(81)-X(82,i)*x(82)-X(83,i)*x(83)-X(84,i)*x(84)-X(85,i)*x(85)-X(86,i)*x(86)-X(87,i)*x(87)-X(88,i)*x(88)-X(89,i)*x(89)-X(90,i)*x(90)-X(91,i)*x(91)-X(92,i)*x(92)-X(93,i)*x(93)-X(94,i)*x(94)-X(95,i)*x(95)-X(96,i)*x(96)-X(97,i)*x(97)-X(98,i)*x(98)-X(99,i)*x(99)-X(100,i)*x(100)-X(101,i)*x(101)-X(102,i)*x(102)-X(103,i)*x(103)-X(104,i)*x(104)-X(105,i)*x(105)-X(106,i)*x(106)-X(107,i)*x(107)-X(108,i)*x(108)-X(109,i)*x(109)-X(110,i)*x(110)-X(111,i)*x(111)-X(112,i)*x(112)-X(113,i)*x(113)-X(114,i)*x(114)-X(115,i)*x(115)-X(116,i)*x(116)-X(117,i)*x(117)-X(118,i)*x(118)-X(119,i)*x(119)-X(120,i)*x(120)-X(121,i)*x(121)-X(122,i)*x(122)-X(123,i)*x(123)-X(124,i)*x(124)-X(125,i)*x(125)-X(126,i)*x(126)-X(127,i)*x(127)-X(128,i)*x(128)-X(129,i)*x(129)-X(130,i)*x(130)-X(131,i)*x(131)-X(132,i)*x(132)-X(133,i)*x(133)-X(134,i)*x(134)-X(135,i)*x(135)-X(136,i)*x(136)-X(137,i)*x(137)-X(138,i)*x(138)-X(139,i)*x(139)-X(140,i)*x(140)-X(141,i)*x(141)-X(142,i)*x(142)-X(143,i)*x(143)-X(144,i)*x(144)-X(145,i)*x(145)-X(146,i)*x(146)-X(147,i)*x(147)-X(148,i)*x(148)-X(149,i)*x(149)-X(150,i)*x(150)-X(151,i)*x(151)-X(152,i)*x(152)-X(153,i)*x(153)-X(154,i)*x(154)-X(155,i)*x(155)-X(156,i)*x(156)-X(157,i)*x(157)-X(158,i)*x(158)-X(159,i)*x(159)-X(160,i)*x(160)-X(161,i)*x(161)-X(162,i)*x(162)-X(163,i)*x(163)-X(164,i)*x(164)-X(165,i)*x(165))*norm(Nvalues(i)))^2+...
            ((u(2,i)-X(1,i)*x(166)-X(2,i)*x(167)-X(3,i)*x(168)-X(4,i)*x(169)-X(5,i)*x(170)-X(6,i)*x(171)-X(7,i)*x(172)-X(8,i)*x(173)-X(9,i)*x(174)-X(10,i)*x(175)-X(11,i)*x(176)-X(12,i)*x(177)-X(13,i)*x(178)-X(14,i)*x(179)-X(15,i)*x(180)-X(16,i)*x(181)-X(17,i)*x(182)-X(18,i)*x(183)-X(19,i)*x(184)-X(20,i)*x(185)-X(21,i)*x(186)-X(22,i)*x(187)-X(23,i)*x(188)-X(24,i)*x(189)-X(25,i)*x(190)-X(26,i)*x(191)-X(27,i)*x(192)-X(28,i)*x(193)-X(29,i)*x(194)-X(30,i)*x(195)-X(31,i)*x(196)-X(32,i)*x(197)-X(33,i)*x(198)-X(34,i)*x(199)-X(35,i)*x(200)-X(36,i)*x(201)-X(37,i)*x(202)-X(38,i)*x(203)-X(39,i)*x(204)-X(40,i)*x(205)-X(41,i)*x(206)-X(42,i)*x(207)-X(43,i)*x(208)-X(44,i)*x(209)-X(45,i)*x(210)-X(46,i)*x(211)-X(47,i)*x(212)-X(48,i)*x(213)-X(49,i)*x(214)-X(50,i)*x(215)-X(51,i)*x(216)-X(52,i)*x(217)-X(53,i)*x(218)-X(54,i)*x(219)-X(55,i)*x(220)-X(56,i)*x(221)-X(57,i)*x(222)-X(58,i)*x(223)-X(59,i)*x(224)-X(60,i)*x(225)-X(61,i)*x(226)-X(62,i)*x(227)-X(63,i)*x(228)-X(64,i)*x(229)-X(65,i)*x(230)-X(66,i)*x(231)-X(67,i)*x(232)-X(68,i)*x(233)-X(69,i)*x(234)-X(70,i)*x(235)-X(71,i)*x(236)-X(72,i)*x(237)-X(73,i)*x(238)-X(74,i)*x(239)-X(75,i)*x(240)-X(76,i)*x(241)-X(77,i)*x(242)-X(78,i)*x(243)-X(79,i)*x(244)-X(80,i)*x(245)-X(81,i)*x(246)-X(82,i)*x(247)-X(83,i)*x(248)-X(84,i)*x(249)-X(85,i)*x(250)-X(86,i)*x(251)-X(87,i)*x(252)-X(88,i)*x(253)-X(89,i)*x(254)-X(90,i)*x(255)-X(91,i)*x(256)-X(92,i)*x(257)-X(93,i)*x(258)-X(94,i)*x(259)-X(95,i)*x(260)-X(96,i)*x(261)-X(97,i)*x(262)-X(98,i)*x(263)-X(99,i)*x(264)-X(100,i)*x(265)-X(101,i)*x(266)-X(102,i)*x(267)-X(103,i)*x(268)-X(104,i)*x(269)-X(105,i)*x(270)-X(106,i)*x(271)-X(107,i)*x(272)-X(108,i)*x(273)-X(109,i)*x(274)-X(110,i)*x(275)-X(111,i)*x(276)-X(112,i)*x(277)-X(113,i)*x(278)-X(114,i)*x(279)-X(115,i)*x(280)-X(116,i)*x(281)-X(117,i)*x(282)-X(118,i)*x(283)-X(119,i)*x(284)-X(120,i)*x(285)-X(121,i)*x(286)-X(122,i)*x(287)-X(123,i)*x(288)-X(124,i)*x(289)-X(125,i)*x(290)-X(126,i)*x(291)-X(127,i)*x(292)-X(128,i)*x(293)-X(129,i)*x(294)-X(130,i)*x(295)-X(131,i)*x(296)-X(132,i)*x(297)-X(133,i)*x(298)-X(134,i)*x(299)-X(135,i)*x(300)-X(136,i)*x(301)-X(137,i)*x(302)-X(138,i)*x(303)-X(139,i)*x(304)-X(140,i)*x(305)-X(141,i)*x(306)-X(142,i)*x(307)-X(143,i)*x(308)-X(144,i)*x(309)-X(145,i)*x(310)-X(146,i)*x(311)-X(147,i)*x(312)-X(148,i)*x(313)-X(149,i)*x(314)-X(150,i)*x(315)-X(151,i)*x(316)-X(152,i)*x(317)-X(153,i)*x(318)-X(154,i)*x(319)-X(155,i)*x(320)-X(156,i)*x(321)-X(157,i)*x(322)-X(158,i)*x(323)-X(159,i)*x(324)-X(160,i)*x(325)-X(161,i)*x(326)-X(162,i)*x(327)-X(163,i)*x(328)-X(164,i)*x(329)-X(165,i)*x(330))*norm(Nvalues(i)))^2+...
            ((u(3,i)-X(1,i)*x(331)-X(2,i)*x(332)-X(3,i)*x(333)-X(4,i)*x(334)-X(5,i)*x(335)-X(6,i)*x(336)-X(7,i)*x(337)-X(8,i)*x(338)-X(9,i)*x(339)-X(10,i)*x(340)-X(11,i)*x(341)-X(12,i)*x(342)-X(13,i)*x(343)-X(14,i)*x(344)-X(15,i)*x(345)-X(16,i)*x(346)-X(17,i)*x(347)-X(18,i)*x(348)-X(19,i)*x(349)-X(20,i)*x(350)-X(21,i)*x(351)-X(22,i)*x(352)-X(23,i)*x(353)-X(24,i)*x(354)-X(25,i)*x(355)-X(26,i)*x(356)-X(27,i)*x(357)-X(28,i)*x(358)-X(29,i)*x(359)-X(30,i)*x(360)-X(31,i)*x(361)-X(32,i)*x(362)-X(33,i)*x(363)-X(34,i)*x(364)-X(35,i)*x(365)-X(36,i)*x(366)-X(37,i)*x(367)-X(38,i)*x(368)-X(39,i)*x(369)-X(40,i)*x(370)-X(41,i)*x(371)-X(42,i)*x(372)-X(43,i)*x(373)-X(44,i)*x(374)-X(45,i)*x(375)-X(46,i)*x(376)-X(47,i)*x(377)-X(48,i)*x(378)-X(49,i)*x(379)-X(50,i)*x(380)-X(51,i)*x(381)-X(52,i)*x(382)-X(53,i)*x(383)-X(54,i)*x(384)-X(55,i)*x(385)-X(56,i)*x(386)-X(57,i)*x(387)-X(58,i)*x(388)-X(59,i)*x(389)-X(60,i)*x(390)-X(61,i)*x(391)-X(62,i)*x(392)-X(63,i)*x(393)-X(64,i)*x(394)-X(65,i)*x(395)-X(66,i)*x(396)-X(67,i)*x(397)-X(68,i)*x(398)-X(69,i)*x(399)-X(70,i)*x(400)-X(71,i)*x(401)-X(72,i)*x(402)-X(73,i)*x(403)-X(74,i)*x(404)-X(75,i)*x(405)-X(76,i)*x(406)-X(77,i)*x(407)-X(78,i)*x(408)-X(79,i)*x(409)-X(80,i)*x(410)-X(81,i)*x(411)-X(82,i)*x(412)-X(83,i)*x(413)-X(84,i)*x(414)-X(85,i)*x(415)-X(86,i)*x(416)-X(87,i)*x(417)-X(88,i)*x(418)-X(89,i)*x(419)-X(90,i)*x(420)-X(91,i)*x(421)-X(92,i)*x(422)-X(93,i)*x(423)-X(94,i)*x(424)-X(95,i)*x(425)-X(96,i)*x(426)-X(97,i)*x(427)-X(98,i)*x(428)-X(99,i)*x(429)-X(100,i)*x(430)-X(101,i)*x(431)-X(102,i)*x(432)-X(103,i)*x(433)-X(104,i)*x(434)-X(105,i)*x(435)-X(106,i)*x(436)-X(107,i)*x(437)-X(108,i)*x(438)-X(109,i)*x(439)-X(110,i)*x(440)-X(111,i)*x(441)-X(112,i)*x(442)-X(113,i)*x(443)-X(114,i)*x(444)-X(115,i)*x(445)-X(116,i)*x(446)-X(117,i)*x(447)-X(118,i)*x(448)-X(119,i)*x(449)-X(120,i)*x(450)-X(121,i)*x(451)-X(122,i)*x(452)-X(123,i)*x(453)-X(124,i)*x(454)-X(125,i)*x(455)-X(126,i)*x(456)-X(127,i)*x(457)-X(128,i)*x(458)-X(129,i)*x(459)-X(130,i)*x(460)-X(131,i)*x(461)-X(132,i)*x(462)-X(133,i)*x(463)-X(134,i)*x(464)-X(135,i)*x(465)-X(136,i)*x(466)-X(137,i)*x(467)-X(138,i)*x(468)-X(139,i)*x(469)-X(140,i)*x(470)-X(141,i)*x(471)-X(142,i)*x(472)-X(143,i)*x(473)-X(144,i)*x(474)-X(145,i)*x(475)-X(146,i)*x(476)-X(147,i)*x(477)-X(148,i)*x(478)-X(149,i)*x(479)-X(150,i)*x(480)-X(151,i)*x(481)-X(152,i)*x(482)-X(153,i)*x(483)-X(154,i)*x(484)-X(155,i)*x(485)-X(156,i)*x(486)-X(157,i)*x(487)-X(158,i)*x(488)-X(159,i)*x(489)-X(160,i)*x(490)-X(161,i)*x(491)-X(162,i)*x(492)-X(163,i)*x(493)-X(164,i)*x(494)-X(165,i)*x(495))*norm(Nvalues(i)))^2;
end

end

end
