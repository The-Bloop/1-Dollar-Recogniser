classdef Dollar1
    methods
        
        function sampled = Sample_Input(obj, points)
            d = diff(points);
            L =cumsum(sqrt(d(:,1).^2 + d(:,2).^2));
            L = [0;L];
            sL = L(end)/63;
            sampled = ones(64,2);
            sampled(1,:) = points(1,:);
            sampled(64,:) = points(end,:);
            sid = [0,0];
            
            sampled(1,:) = points(1,:);
            sampled(end,:) = points(end,:);
            
            for i=1:62
                id = find(L>(i*sL)); 
                sid(1) = id(1) - 1;
                sid(2) = id(1);
                dl = L(sid(2)) - i*sL;
                l = L(sid(2)) - L(sid(1));
                tl = l-dl;
                a = (tl/l);
                sampled(i+1,:) = a.*d(sid(1),:) + points(sid(1),:);
            end
            x = diff(sampled);
            difference = sqrt(x(:,1).^2 + x(:,2).^2);
            
        end
        
        function [centroid,transformed] = Create_Template(obj, sampled)
           centroid = mean(sampled);
           start_centroid = centroid - sampled(1,:);
           theta = atan2(start_centroid(2),start_centroid(1));
           T1 = eye(3); 
           T1(1,3) = -centroid(1);     
           T1(2,3)= -centroid(2);
           T2 = [cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
           z = ones(64,1);
           
           sampled = [sampled z];
           size(sampled);
           transformed = (T2 * T1 * sampled')';
           
           bb_max = max(transformed,[],1);
           bb_min = min(transformed,[],1);
           range = bb_max - bb_min;
           sx = 1 / range(1);
           sy = 1 / range(2);
           T3 = eye(3);
           T3(1,1) = sx;
           T3(2,2) = sy;
           transformed = (T3 * (transformed'))';
           
           bb_min_new = min(transformed,[],1);
           T4 = eye(3);
           T4(1,3) = -bb_min_new(1);
           T4(2,3) = -bb_min_new(2);
           transformed = (T4 * (transformed'))';
           
           new_centroid = (T4 * T3 * T2 * T1 * [centroid 1]')';
           
        end
        
        function index = Get_Match(obj,C,T)
            [m,n,o] = size(T);
            score = zeros(o,1);
            for i=1:o
                score(i) = obj.Match_Template(C,T(:,:,i));
            end
            max_score = max(score);
            index = find(score == max_score);
        end
        
        function flag = Match_Template(obj,C,T)
            
            dif = C - T;
            d = sum((sqrt(dif(:,1).^2 + dif(:,2).^2)))/32;
            score = 1 - ((2*d)/sqrt(2));
            
            dtheta = pi/180;
            
            R1 = [cos(dtheta),sin(dtheta),0;-sin(dtheta),cos(dtheta),0;0,0,1];
            R2 = [cos(-dtheta),sin(-dtheta),0;-sin(-dtheta),cos(-dtheta),0;0,0,1];
            
            T1 = (R1 * T')';
            T2 = (R2 * T')';
            
            dif1 = C - T1;
            d1 = sum((sqrt(dif1(:,1).^2 + dif1(:,2).^2)))/32;
            score1 = 1 - ((2*d1)/sqrt(2));
            
            dif2 = C - T2;
            d2 = sum((sqrt(dif2(:,1).^2 + dif2(:,2).^2)))/32;
            score2 = 1 - ((2*d2)/sqrt(2));
            
            theta = 0;
            score_max = score;
            if(score1 > score)
                dtheta = pi/180;
                theta = dtheta * 2;
                score_max = score1;
            elseif(score2 > score)
                dtheta = -pi/180;
                theta = -dtheta * 2;
                score_max = score2;
            end

            for i = 1:10
                Rn = [cos(dtheta),sin(dtheta),0;-sin(dtheta),cos(dtheta),0;0,0,1];
                Tn = (Rn * T')';
                difn = C - Tn;
                dn = sum((sqrt(difn(:,1).^2 + difn(:,2).^2)))/32;
                scoren = 1 - ((2*dn)/sqrt(2));
                if(scoren > score_max)
                    score_max = scoren;
                end
                theta = theta + dtheta;
            end
            flag = score_max;
        end
        
    end
end