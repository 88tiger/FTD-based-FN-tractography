function Trace = ftrack(A,start_point,step_size,mask,RGVPregion,startpoint)

%% Parameter
detaT = 0.1;
stop = 0;
reach = 0;
p0 = start_point;
Trace = zeros(0,3); 
stopTrace = zeros(0,3);
oldvector = [0 0 0];
checkintp0 = [0 0 0];

%% Tracking
while(1)
    intp0 = ceil(p0);
    
    % Maximum length
    if size(Trace,1)>100/step_size
        break
    end
    % Out of size
    if any(intp0(:,3)<1) || any(intp0(:,3)>size(mask,3)) || ...
            any(intp0(:,2)<1) || any(intp0(:,2)>size(mask,2)) || ...
            any(intp0(:,1)<1) || any(intp0(:,1)>size(mask,1))
%         Trace = [Trace;stopTrace;p0]; stopTrace = [];
        break;
    end
    
    % Reach end point
    if startpoint(intp0(1),intp0(2),intp0(3))==1 && reach==0
        reach = 1;
    elseif startpoint(intp0(1),intp0(2),intp0(3))==0 && reach==1 
        break;
    end
    
    
    % Out of mask
    if RGVPregion(intp0(1),intp0(2),intp0(3))==0 && any(checkintp0~=intp0)
        stop = stop+1;
    elseif RGVPregion(intp0(1),intp0(2),intp0(3))==1 && stop>0
        stop = 0;
    end
    if stop>=1
        break
    end
    checkintp0 = intp0;
    
%     if mask(intp0(1),intp0(2),intp0(3)) == 0
% %         Trace = [Trace;stopTrace;p0]; stopTrace = [];
%         break;
%     end
    
    % Cal A
    if size(A,2) == 4
        f = @(p) (A*[p(1), p(2), p(3), 1]')';
    elseif size(A,2) == 6
        f = @(p) (A*[p(1)^2, p(2)^2, p(1)*p(2), p(1), p(2), 1]')';
    elseif size(A,2) == 10 && ~FiberCup
        f = @(p) (A*[p(1)^2, p(2)^2, p(3)^2, p(1)*p(2), p(1)*p(3), p(2)*p(3), p(1), p(2), p(3), 1]')';
    elseif size(A,2) == 10 && FiberCup
        f = @(p) (A*[p(1)^3, p(2)^3, (p(1)^2)*p(2), p(1)*(p(2)^2), p(1)^2, p(2)^2, p(1)*p(2), p(1), p(2), 1]')';
    elseif size(A,2) == 15
        f = @(p) (A*[p(1)^4, p(2)^4, (p(1)^3)*p(2), p(1)*(p(2)^3), (p(1)^2)*(p(2)^2),p(1)^3, p(2)^3, (p(1)^2)*p(2), p(1)*(p(2)^2), p(1)^2, p(2)^2, ...
                    p(1)*p(2), p(1), p(2), 1]')';  
    elseif size(A,2) == 20
        f = @(p) (A*[p(1)^3, p(2)^3, p(3)^3, (p(1)^2)*p(2), (p(1)^2)*p(3), (p(2)^2)*p(1), (p(2)^2)*p(3), (p(3)^2)*p(1), (p(3)^2)*p(2), p(1)*p(2)*p(3), ...
                    p(1)^2, p(2)^2, p(3)^2, p(1)*p(2), p(1)*p(3), p(2)*p(3), p(1), p(2), p(3), 1]')';
    elseif size(A,2) == 35
        f = @(p) (A*[p(1)^4, p(2)^4, p(3)^4, (p(1)^3)*p(2), (p(1)^3)*p(3), (p(2)^3)*p(1), (p(2)^3)*p(3), (p(3)^3)*p(1), (p(3)^3)*p(2), (p(1)^2)*(p(2)^2), (p(1)^2)*(p(3)^2), (p(2)^2)*(p(3)^2), (p(1)^2)*p(2)*p(3), (p(2)^2)*p(1)*p(3), (p(3)^2)*p(1)*p(2), ...
                    p(1)^3, p(2)^3, p(3)^3, (p(1)^2)*p(2), (p(1)^2)*p(3), (p(2)^2)*p(1), (p(2)^2)*p(3), (p(3)^2)*p(1), (p(3)^2)*p(2), p(1)*p(2)*p(3), ...
                    p(1)^2, p(2)^2, p(3)^2, p(1)*p(2), p(1)*p(3), p(2)*p(3), ...
                    p(1), p(2), p(3), 1]')';
    elseif size(A,2) == 56
        f = @(p) (A*[p(1)^5, p(2)^5, p(3)^5, (p(1)^4)*p(2), (p(1)^4)*p(3), (p(2)^4)*p(1), (p(2)^4)*p(3), (p(3)^4)*p(1), (p(3)^4)*p(2), (p(1)^3)*(p(2)^2), (p(1)^3)*(p(3)^2), (p(2)^3)*(p(1)^2), (p(2)^3)*(p(3)^2), (p(3)^3)*(p(1)^2), (p(3)^3)*(p(2)^2), (p(1)^3)*p(2)*p(3), (p(2)^3)*p(1)*p(3), (p(3)^3)*p(1)*p(2), (p(1)^2)*(p(2)^2)*p(3), (p(1)^2)*p(2)*(p(3)^2), p(1)*(p(2)^2)*(p(3)^2), ...
                    p(1)^4, p(2)^4, p(3)^4, (p(1)^3)*p(2), (p(1)^3)*p(3), (p(2)^3)*p(1), (p(2)^3)*p(3), (p(3)^3)*p(1), (p(3)^3)*p(2), (p(1)^2)*(p(2)^2), (p(1)^2)*(p(3)^2), (p(2)^2)*(p(3)^2), (p(1)^2)*p(2)*p(3), (p(2)^2)*p(1)*p(3), (p(3)^2)*p(1)*p(2), ...
                    p(1)^3, p(2)^3, p(3)^3, (p(1)^2)*p(2), (p(1)^2)*p(3), (p(2)^2)*p(1), (p(2)^2)*p(3), (p(3)^2)*p(1), (p(3)^2)*p(2), p(1)*p(2)*p(3), ...
                    p(1)^2, p(2)^2, p(3)^2, p(1)*p(2), p(1)*p(3), p(2)*p(3), ...
                    p(1), p(2), p(3), 1]')';
    elseif size(A,2) == 84
        f = @(p) (A*[p(1), p(1)^2, p(1)^3, p(1)^4, p(1)^5, p(1)^6, p(2), p(1)*p(2), p(1)^2*p(2), p(1)^3*p(2), p(1)^4*p(2), p(1)^5*p(2), ...
                    p(2)^2, p(1)*p(2)^2, p(1)^2*p(2)^2, p(1)^3*p(2)^2, p(1)^4*p(2)^2, p(2)^3, p(1)*p(2)^3, p(1)^2*p(2)^3, p(1)^3*p(2)^3, ...
                    p(2)^4, p(1)*p(2)^4, p(1)^2*p(2)^4, p(2)^5, p(1)*p(2)^5, p(2)^6, p(3), p(1)*p(3), p(1)^2*p(3), p(1)^3*p(3), p(1)^4*p(3), p(1)^5*p(3), ...
                    p(2)*p(3), p(1)*p(2)*p(3), p(1)^2*p(2)*p(3), p(1)^3*p(2)*p(3), p(1)^4*p(2)*p(3), p(2)^2*p(3), p(1)*p(2)^2*p(3), p(1)^2*p(2)^2*p(3), p(1)^3*p(2)^2*p(3), ...
                    p(2)^3*p(3), p(1)*p(2)^3*p(3), p(1)^2*p(2)^3*p(3), p(2)^4*p(3), p(1)*p(2)^4*p(3), p(2)^5*p(3), p(3)^2, p(1)*p(3)^2, p(1)^2*p(3)^2, p(1)^3*p(3)^2, ...
                    p(1)^4*p(3)^2, p(2)*p(3)^2, p(1)*p(2)*p(3)^2, p(1)^2*p(2)*p(3)^2, p(1)^3*p(2)*p(3)^2, p(2)^2*p(3)^2, p(1)*p(2)^2*p(3)^2, p(1)^2*p(2)^2*p(3)^2, ...
                    p(2)^3*p(3)^2, p(1)*p(2)^3*p(3)^2, p(2)^4*p(3)^2, p(3)^3, p(1)*p(3)^3, p(1)^2*p(3)^3, p(1)^3*p(3)^3, p(2)*p(3)^3, p(1)*p(2)*p(3)^3, p(1)^2*p(2)*p(3)^3, ...
                    p(2)^2*p(3)^3, p(1)*p(2)^2*p(3)^3, p(2)^3*p(3)^3,p(3)^4, p(1)*p(3)^4, p(1)^2*p(3)^4, p(2)*p(3)^4, p(1)*p(2)*p(3)^4, p(2)^2*p(3)^4, p(3)^5, p(1)*p(3)^5, p(2)*p(3)^5, p(3)^6, 1]')';
    elseif size(A,2) == 120
        f = @(p) (A*[p(1), p(1)^2, p(1)^3, p(1)^4, p(1)^5, p(1)^6, p(1)^7, p(2), p(1)*p(2), p(1)^2*p(2), p(1)^3*p(2), p(1)^4*p(2), p(1)^5*p(2), p(1)^6*p(2), p(2)^2, p(1)*p(2)^2, ...
                    p(1)^2*p(2)^2, p(1)^3*p(2)^2, p(1)^4*p(2)^2, p(1)^5*p(2)^2, p(2)^3, p(1)*p(2)^3, p(1)^2*p(2)^3, p(1)^3*p(2)^3, p(1)^4*p(2)^3, p(2)^4, p(1)*p(2)^4, p(1)^2*p(2)^4, ...
                    p(1)^3*p(2)^4, p(2)^5, p(1)*p(2)^5, p(1)^2*p(2)^5, p(2)^6, p(1)*p(2)^6, p(2)^7, p(3), p(1)*p(3), p(1)^2*p(3), p(1)^3*p(3), p(1)^4*p(3), p(1)^5*p(3), p(1)^6*p(3), p(2)*p(3), ...
                    p(1)*p(2)*p(3), p(1)^2*p(2)*p(3), p(1)^3*p(2)*p(3), p(1)^4*p(2)*p(3), p(1)^5*p(2)*p(3), p(2)^2*p(3), p(1)*p(2)^2*p(3), p(1)^2*p(2)^2*p(3), p(1)^3*p(2)^2*p(3), p(1)^4*p(2)^2*p(3), ...
                    p(2)^3*p(3), p(1)*p(2)^3*p(3), p(1)^2*p(2)^3*p(3), p(1)^3*p(2)^3*p(3), p(2)^4*p(3), p(1)*p(2)^4*p(3), p(1)^2*p(2)^4*p(3), p(2)^5*p(3), p(1)*p(2)^5*p(3), p(2)^6*p(3), p(3)^2, ...
                    p(1)*p(3)^2, p(1)^2*p(3)^2, p(1)^3*p(3)^2, p(1)^4*p(3)^2, p(1)^5*p(3)^2, p(2)*p(3)^2, p(1)*p(2)*p(3)^2, p(1)^2*p(2)*p(3)^2, p(1)^3*p(2)*p(3)^2, p(1)^4*p(2)*p(3)^2, ...
                    p(2)^2*p(3)^2, p(1)*p(2)^2*p(3)^2, p(1)^2*p(2)^2*p(3)^2, p(1)^3*p(2)^2*p(3)^2, p(2)^3*p(3)^2, p(1)*p(2)^3*p(3)^2, p(1)^2*p(2)^3*p(3)^2, p(2)^4*p(3)^2, p(1)*p(2)^4*p(3)^2, ...
                    p(2)^5*p(3)^2, p(3)^3, p(1)*p(3)^3, p(1)^2*p(3)^3, p(1)^3*p(3)^3, p(1)^4*p(3)^3, p(2)*p(3)^3, p(1)*p(2)*p(3)^3, p(1)^2*p(2)*p(3)^3, p(1)^3*p(2)*p(3)^3, p(2)^2*p(3)^3, ...
                    p(1)*p(2)^2*p(3)^3, p(1)^2*p(2)^2*p(3)^3, p(2)^3*p(3)^3, p(1)*p(2)^3*p(3)^3, p(2)^4*p(3)^3, p(3)^4, p(1)*p(3)^4, p(1)^2*p(3)^4, p(1)^3*p(3)^4, p(2)*p(3)^4, p(1)*p(2)*p(3)^4, ...
                    p(1)^2*p(2)*p(3)^4, p(2)^2*p(3)^4, p(1)*p(2)^2*p(3)^4, p(2)^3*p(3)^4, p(3)^5, p(1)*p(3)^5, p(1)^2*p(3)^5, p(2)*p(3)^5, p(1)*p(2)*p(3)^5, p(2)^2*p(3)^5, p(3)^6, p(1)*p(3)^6, ...
                    p(2)*p(3)^6, p(3)^7, 1]')';
    elseif size(A,2) == 165
        f = @(p) (A*[p(1), p(1)^2, p(1)^3, p(1)^4, p(1)^5, p(1)^6, p(1)^7, p(1)^8, p(2), p(1)*p(2), p(1)^2*p(2), p(1)^3*p(2), p(1)^4*p(2), p(1)^5*p(2), p(1)^6*p(2), p(1)^7*p(2), ...
                    p(2)^2, p(1)*p(2)^2, p(1)^2*p(2)^2, p(1)^3*p(2)^2, p(1)^4*p(2)^2, p(1)^5*p(2)^2, p(1)^6*p(2)^2, p(2)^3, p(1)*p(2)^3, p(1)^2*p(2)^3, p(1)^3*p(2)^3, p(1)^4*p(2)^3, ...
                    p(1)^5*p(2)^3, p(2)^4, p(1)*p(2)^4, p(1)^2*p(2)^4, p(1)^3*p(2)^4, p(1)^4*p(2)^4, p(2)^5, p(1)*p(2)^5, p(1)^2*p(2)^5, p(1)^3*p(2)^5, p(2)^6, p(1)*p(2)^6, p(1)^2*p(2)^6, ...
                    p(2)^7, p(1)*p(2)^7, p(2)^8, p(3), p(1)*p(3), p(1)^2*p(3), p(1)^3*p(3), p(1)^4*p(3), p(1)^5*p(3), p(1)^6*p(3), p(1)^7*p(3), p(2)*p(3), p(1)*p(2)*p(3), p(1)^2*p(2)*p(3), ...
                    p(1)^3*p(2)*p(3), p(1)^4*p(2)*p(3), p(1)^5*p(2)*p(3), p(1)^6*p(2)*p(3), p(2)^2*p(3), p(1)*p(2)^2*p(3), p(1)^2*p(2)^2*p(3), p(1)^3*p(2)^2*p(3), p(1)^4*p(2)^2*p(3), p(1)^5*p(2)^2*p(3), ...
                    p(2)^3*p(3), p(1)*p(2)^3*p(3), p(1)^2*p(2)^3*p(3), p(1)^3*p(2)^3*p(3), p(1)^4*p(2)^3*p(3), p(2)^4*p(3), p(1)*p(2)^4*p(3), p(1)^2*p(2)^4*p(3), p(1)^3*p(2)^4*p(3), p(2)^5*p(3), ...
                    p(1)*p(2)^5*p(3), p(1)^2*p(2)^5*p(3), p(2)^6*p(3), p(1)*p(2)^6*p(3), p(2)^7*p(3), p(3)^2, p(1)*p(3)^2, p(1)^2*p(3)^2, p(1)^3*p(3)^2, p(1)^4*p(3)^2, p(1)^5*p(3)^2, p(1)^6*p(3)^2, ...
                    p(2)*p(3)^2, p(1)*p(2)*p(3)^2, p(1)^2*p(2)*p(3)^2, p(1)^3*p(2)*p(3)^2, p(1)^4*p(2)*p(3)^2, p(1)^5*p(2)*p(3)^2, p(2)^2*p(3)^2, p(1)*p(2)^2*p(3)^2, p(1)^2*p(2)^2*p(3)^2, ...
                    p(1)^3*p(2)^2*p(3)^2, p(1)^4*p(2)^2*p(3)^2, p(2)^3*p(3)^2, p(1)*p(2)^3*p(3)^2, p(1)^2*p(2)^3*p(3)^2, p(1)^3*p(2)^3*p(3)^2, p(2)^4*p(3)^2, p(1)*p(2)^4*p(3)^2, ...
                    p(1)^2*p(2)^4*p(3)^2, p(2)^5*p(3)^2, p(1)*p(2)^5*p(3)^2, p(2)^6*p(3)^2, p(3)^3, p(1)*p(3)^3, p(1)^2*p(3)^3, p(1)^3*p(3)^3, p(1)^4*p(3)^3, p(1)^5*p(3)^3, p(2)*p(3)^3, ...
                    p(1)*p(2)*p(3)^3, p(1)^2*p(2)*p(3)^3, p(1)^3*p(2)*p(3)^3, p(1)^4*p(2)*p(3)^3, p(2)^2*p(3)^3, p(1)*p(2)^2*p(3)^3, p(1)^2*p(2)^2*p(3)^3, p(1)^3*p(2)^2*p(3)^3, p(2)^3*p(3)^3, ...
                    p(1)*p(2)^3*p(3)^3, p(1)^2*p(2)^3*p(3)^3, p(2)^4*p(3)^3, p(1)*p(2)^4*p(3)^3, p(2)^5*p(3)^3, p(3)^4, p(1)*p(3)^4, p(1)^2*p(3)^4, p(1)^3*p(3)^4, p(1)^4*p(3)^4, p(2)*p(3)^4, ...
                    p(1)*p(2)*p(3)^4, p(1)^2*p(2)*p(3)^4, p(1)^3*p(2)*p(3)^4, p(2)^2*p(3)^4, p(1)*p(2)^2*p(3)^4, p(1)^2*p(2)^2*p(3)^4, p(2)^3*p(3)^4, p(1)*p(2)^3*p(3)^4, p(2)^4*p(3)^4, p(3)^5, ...
                    p(1)*p(3)^5, p(1)^2*p(3)^5, p(1)^3*p(3)^5, p(2)*p(3)^5, p(1)*p(2)*p(3)^5, p(1)^2*p(2)*p(3)^5, p(2)^2*p(3)^5, p(1)*p(2)^2*p(3)^5, p(2)^3*p(3)^5, p(3)^6, p(1)*p(3)^6, ...
                    p(1)^2*p(3)^6, p(2)*p(3)^6, p(1)*p(2)*p(3)^6, p(2)^2*p(3)^6, p(3)^7, p(1)*p(3)^7, p(2)*p(3)^7, p(3)^8, 1]')';
    elseif size(A,2) == 286
        f = @(p) (A*[p(1), p(1)^2, p(1)^3, p(1)^4, p(1)^5, p(1)^6, p(1)^7, p(1)^8, p(1)^9, p(1)^10, p(2), p(1)*p(2), p(1)^2*p(2), p(1)^3*p(2), ...
                    p(1)^4*p(2), p(1)^5*p(2), p(1)^6*p(2), p(1)^7*p(2), p(1)^8*p(2), p(1)^9*p(2), p(2)^2, p(1)*p(2)^2, p(1)^2*p(2)^2, p(1)^3*p(2)^2, ...
                    p(1)^4*p(2)^2, p(1)^5*p(2)^2, p(1)^6*p(2)^2, p(1)^7*p(2)^2, p(1)^8*p(2)^2, p(2)^3, p(1)*p(2)^3, p(1)^2*p(2)^3, p(1)^3*p(2)^3, ...
                    p(1)^4*p(2)^3, p(1)^5*p(2)^3, p(1)^6*p(2)^3, p(1)^7*p(2)^3, p(2)^4, p(1)*p(2)^4, p(1)^2*p(2)^4, p(1)^3*p(2)^4, p(1)^4*p(2)^4, ...
                    p(1)^5*p(2)^4, p(1)^6*p(2)^4, p(2)^5, p(1)*p(2)^5, p(1)^2*p(2)^5, p(1)^3*p(2)^5, p(1)^4*p(2)^5, p(1)^5*p(2)^5, p(2)^6, p(1)*p(2)^6, ...
                    p(1)^2*p(2)^6, p(1)^3*p(2)^6, p(1)^4*p(2)^6, p(2)^7, p(1)*p(2)^7, p(1)^2*p(2)^7, p(1)^3*p(2)^7, p(2)^8, p(1)*p(2)^8, p(1)^2*p(2)^8, ...
                    p(2)^9, p(1)*p(2)^9, p(2)^10, p(3), p(1)*p(3), p(1)^2*p(3), p(1)^3*p(3), p(1)^4*p(3), p(1)^5*p(3), p(1)^6*p(3), p(1)^7*p(3), ...
                    p(1)^8*p(3), p(1)^9*p(3), p(2)*p(3), p(1)*p(2)*p(3), p(1)^2*p(2)*p(3), p(1)^3*p(2)*p(3), p(1)^4*p(2)*p(3), p(1)^5*p(2)*p(3), ...
                    p(1)^6*p(2)*p(3), p(1)^7*p(2)*p(3), p(1)^8*p(2)*p(3), p(2)^2*p(3), p(1)*p(2)^2*p(3), p(1)^2*p(2)^2*p(3), p(1)^3*p(2)^2*p(3), ...
                    p(1)^4*p(2)^2*p(3), p(1)^5*p(2)^2*p(3), p(1)^6*p(2)^2*p(3), p(1)^7*p(2)^2*p(3), p(2)^3*p(3), p(1)*p(2)^3*p(3), p(1)^2*p(2)^3*p(3), ...
                    p(1)^3*p(2)^3*p(3), p(1)^4*p(2)^3*p(3), p(1)^5*p(2)^3*p(3), p(1)^6*p(2)^3*p(3), p(2)^4*p(3), p(1)*p(2)^4*p(3), p(1)^2*p(2)^4*p(3), ...
                    p(1)^3*p(2)^4*p(3), p(1)^4*p(2)^4*p(3), p(1)^5*p(2)^4*p(3), p(2)^5*p(3), p(1)*p(2)^5*p(3), p(1)^2*p(2)^5*p(3), p(1)^3*p(2)^5*p(3), ...
                    p(1)^4*p(2)^5*p(3), p(2)^6*p(3), p(1)*p(2)^6*p(3), p(1)^2*p(2)^6*p(3), p(1)^3*p(2)^6*p(3), p(2)^7*p(3), p(1)*p(2)^7*p(3), ...
                    p(1)^2*p(2)^7*p(3), p(2)^8*p(3), p(1)*p(2)^8*p(3), p(2)^9*p(3), p(3)^2, p(1)*p(3)^2, p(1)^2*p(3)^2, p(1)^3*p(3)^2, p(1)^4*p(3)^2, ...
                    p(1)^5*p(3)^2, p(1)^6*p(3)^2, p(1)^7*p(3)^2, p(1)^8*p(3)^2, p(2)*p(3)^2, p(1)*p(2)*p(3)^2, p(1)^2*p(2)*p(3)^2, p(1)^3*p(2)*p(3)^2, ...
                    p(1)^4*p(2)*p(3)^2, p(1)^5*p(2)*p(3)^2, p(1)^6*p(2)*p(3)^2, p(1)^7*p(2)*p(3)^2, p(2)^2*p(3)^2, p(1)*p(2)^2*p(3)^2, ...
                    p(1)^2*p(2)^2*p(3)^2, p(1)^3*p(2)^2*p(3)^2, p(1)^4*p(2)^2*p(3)^2, p(1)^5*p(2)^2*p(3)^2, p(1)^6*p(2)^2*p(3)^2, p(2)^3*p(3)^2, ...
                    p(1)*p(2)^3*p(3)^2, p(1)^2*p(2)^3*p(3)^2, p(1)^3*p(2)^3*p(3)^2, p(1)^4*p(2)^3*p(3)^2, p(1)^5*p(2)^3*p(3)^2, p(2)^4*p(3)^2, ...
                    p(1)*p(2)^4*p(3)^2, p(1)^2*p(2)^4*p(3)^2, p(1)^3*p(2)^4*p(3)^2, p(1)^4*p(2)^4*p(3)^2, p(2)^5*p(3)^2, p(1)*p(2)^5*p(3)^2, ...
                    p(1)^2*p(2)^5*p(3)^2, p(1)^3*p(2)^5*p(3)^2, p(2)^6*p(3)^2, p(1)*p(2)^6*p(3)^2, p(1)^2*p(2)^6*p(3)^2, p(2)^7*p(3)^2, ...
                    p(1)*p(2)^7*p(3)^2, p(2)^8*p(3)^2, p(3)^3, p(1)*p(3)^3, p(1)^2*p(3)^3, p(1)^3*p(3)^3, p(1)^4*p(3)^3, p(1)^5*p(3)^3, ...
                    p(1)^6*p(3)^3, p(1)^7*p(3)^3, p(2)*p(3)^3, p(1)*p(2)*p(3)^3, p(1)^2*p(2)*p(3)^3, p(1)^3*p(2)*p(3)^3, p(1)^4*p(2)*p(3)^3, ...
                    p(1)^5*p(2)*p(3)^3, p(1)^6*p(2)*p(3)^3, p(2)^2*p(3)^3, p(1)*p(2)^2*p(3)^3, p(1)^2*p(2)^2*p(3)^3, p(1)^3*p(2)^2*p(3)^3, ...
                    p(1)^4*p(2)^2*p(3)^3, p(1)^5*p(2)^2*p(3)^3, p(2)^3*p(3)^3, p(1)*p(2)^3*p(3)^3, p(1)^2*p(2)^3*p(3)^3, p(1)^3*p(2)^3*p(3)^3, ...
                    p(1)^4*p(2)^3*p(3)^3, p(2)^4*p(3)^3, p(1)*p(2)^4*p(3)^3, p(1)^2*p(2)^4*p(3)^3, p(1)^3*p(2)^4*p(3)^3, p(2)^5*p(3)^3, ...
                    p(1)*p(2)^5*p(3)^3, p(1)^2*p(2)^5*p(3)^3, p(2)^6*p(3)^3, p(1)*p(2)^6*p(3)^3, p(2)^7*p(3)^3, p(3)^4, p(1)*p(3)^4, p(1)^2*p(3)^4, ...
                    p(1)^3*p(3)^4, p(1)^4*p(3)^4, p(1)^5*p(3)^4, p(1)^6*p(3)^4, p(2)*p(3)^4, p(1)*p(2)*p(3)^4, p(1)^2*p(2)*p(3)^4, p(1)^3*p(2)*p(3)^4, ...
                    p(1)^4*p(2)*p(3)^4, p(1)^5*p(2)*p(3)^4, p(2)^2*p(3)^4, p(1)*p(2)^2*p(3)^4, p(1)^2*p(2)^2*p(3)^4, p(1)^3*p(2)^2*p(3)^4, ...
                    p(1)^4*p(2)^2*p(3)^4, p(2)^3*p(3)^4, p(1)*p(2)^3*p(3)^4, p(1)^2*p(2)^3*p(3)^4, p(1)^3*p(2)^3*p(3)^4, p(2)^4*p(3)^4, ...
                    p(1)*p(2)^4*p(3)^4, p(1)^2*p(2)^4*p(3)^4, p(2)^5*p(3)^4, p(1)*p(2)^5*p(3)^4, p(2)^6*p(3)^4, p(3)^5, p(1)*p(3)^5, p(1)^2*p(3)^5, ...
                    p(1)^3*p(3)^5, p(1)^4*p(3)^5, p(1)^5*p(3)^5, p(2)*p(3)^5, p(1)*p(2)*p(3)^5, p(1)^2*p(2)*p(3)^5, p(1)^3*p(2)*p(3)^5, ...
                    p(1)^4*p(2)*p(3)^5, p(2)^2*p(3)^5, p(1)*p(2)^2*p(3)^5, p(1)^2*p(2)^2*p(3)^5, p(1)^3*p(2)^2*p(3)^5, p(2)^3*p(3)^5, ...
                    p(1)*p(2)^3*p(3)^5, p(1)^2*p(2)^3*p(3)^5, p(2)^4*p(3)^5, p(1)*p(2)^4*p(3)^5, p(2)^5*p(3)^5, p(3)^6, p(1)*p(3)^6, p(1)^2*p(3)^6, ...
                    p(1)^3*p(3)^6, p(1)^4*p(3)^6, p(2)*p(3)^6, p(1)*p(2)*p(3)^6, p(1)^2*p(2)*p(3)^6, p(1)^3*p(2)*p(3)^6, p(2)^2*p(3)^6, ...
                    p(1)*p(2)^2*p(3)^6, p(1)^2*p(2)^2*p(3)^6, p(2)^3*p(3)^6, p(1)*p(2)^3*p(3)^6, p(2)^4*p(3)^6, p(3)^7, p(1)*p(3)^7, ...
                    p(1)^2*p(3)^7, p(1)^3*p(3)^7, p(2)*p(3)^7, p(1)*p(2)*p(3)^7, p(1)^2*p(2)*p(3)^7, p(2)^2*p(3)^7, p(1)*p(2)^2*p(3)^7, ...
                    p(2)^3*p(3)^7, p(3)^8, p(1)*p(3)^8, p(1)^2*p(3)^8, p(2)*p(3)^8, p(1)*p(2)*p(3)^8, p(2)^2*p(3)^8, p(3)^9, p(1)*p(3)^9, p(2)*p(3)^9, ...
                    p(3)^10, 1]')';
    end
    
    % Disturbed flow
    disturbed_vector = runge_kutta_4_vector(f, p0, detaT);

    p0 = p0 + (step_size*disturbed_vector(:)/norm(disturbed_vector))';
    
    % Angle thr
    if ~all(oldvector==0)
%         (oldvector/norm(oldvector))*(disturbed_vector(:)/norm(disturbed_vector))
        if (oldvector/norm(oldvector))*(disturbed_vector(:)/norm(disturbed_vector))<0
%             Trace = [Trace;stopTrace;p0]; stopTrace = [];
%             stopreason{fibernum+1, num} = 'Angle';
            break;
        else
            Trace = [Trace;p0];
        end
    end
    oldvector = disturbed_vector;
end

end