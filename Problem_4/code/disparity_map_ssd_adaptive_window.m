%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates disparity map using SSD for defined window size.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function disparity_map_ssd_adaptive_window(I_left, I_right, window_size)

% Add boundary to apply SSD on the edge pixels
I_left_mod = zeros(window_size-1+size(I_left,1),window_size-1+size(I_left,2));
I_left_mod(ceil(window_size/2):end-floor(window_size/2),ceil(window_size/2):end-floor(window_size/2)) = I_left;
I_right_mod = zeros(window_size-1+size(I_right,1),window_size-1+size(I_right,2));
I_right_mod(ceil(window_size/2):end-floor(window_size/2),ceil(window_size/2):end-floor(window_size/2)) = I_right;

% Compute disparity using SSD
D = zeros(size(I_left));
for i=ceil(window_size/2):floor(window_size/2)+size(I_left,1)
    for j=ceil(window_size/2):floor(window_size/2)+size(I_left,2)
        min_diff = Inf;
        % Create the left window
        left_window = I_left_mod(i-floor(window_size/2):i+floor(window_size/2),j-floor(window_size/2):j+floor(window_size/2));
        for k=ceil(window_size/2):floor(window_size/2)+size(I_left,2)
            % Create the right window
            right_window = I_right_mod(i-floor(window_size/2):i+floor(window_size/2),k-floor(window_size/2):k+floor(window_size/2));
            % Compute SSD
            ssd = sum(sum((left_window - right_window).^2));
            if ssd<min_diff
                min_diff = ssd;
                % If SSD is minimum, store disparity
                D(i-floor(window_size/2),j-floor(window_size/2)) = k - j;
            end
        end
    end
end

% Apply adaptive window
actions = [0 0; 1 0; 0 1; -1 0; 0 -1];
D_improv = D;
for i=ceil(window_size/2):floor(window_size/2)+size(I_left,1)
    for j=ceil(window_size/2):floor(window_size/2)+size(I_left,2)
        % Compute central disparity
        D_central = D(i-floor(window_size/2),j-floor(window_size/2));
        % Define initial window size
        start_window = [3 3];
        start_lower_bound = [1 1];
        start_upper_bound = [2 2];
        % Initiate min_var
        min_var = [];
        % Initiate motion allowed in all directions
        motion_allowed = [1 1 1 1 1];
        while any(motion_allowed(2:5))&&(max(start_upper_bound)<=floor(window_size/2))&&(max(start_lower_bound)<=floor(window_size/2))
            var = Inf*ones(1,4);
            for step=1:5
                if motion_allowed(step)
                    % Get the current window
                    current_window = start_window + actions(step,:);
                    if step==1
                        lower_bound = start_lower_bound;
                        upper_bound = start_upper_bound;
                    elseif (step==2)||(step==3)
                        lower_bound = start_lower_bound;
                        upper_bound = start_upper_bound + actions(step,:);
                    elseif (step==4)||(step==5)
                        lower_bound = start_lower_bound - actions(step,:);
                        upper_bound = start_upper_bound;
                    end
                    % Create the right window
                    right_window = I_right_mod(i-lower_bound(1):i+upper_bound(1)-1,j-lower_bound(2)+D_central:j+upper_bound(2)-1+D_central);
                    % Compute alphas
                    alpha_d = 0;
                    alpha_f = 0;
                    for k=1:current_window(1)
                        for l=1:current_window(2)
                            if l == 1
                                alpha_f = alpha_f + (right_window(k,2) - right_window(k,1))^2;
                            elseif l == current_window(2)
                                alpha_f = alpha_f + (right_window(k,l) - right_window(k,l-1))^2;
                            else
                                alpha_f = alpha_f + ((right_window(k,l+1) - right_window(k,l-1))/2)^2;
                            end
                            try
                                alpha_d = alpha_d + ((D(i-floor(window_size/2)+k-lower_bound(1)-1,j-floor(window_size/2)+l-lower_bound(2)-1) - D_central)^2)/sqrt((i-floor(window_size/2)+k-lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-lower_bound(2)-1)^2);
                            catch
                                continue;
                            end
                        end
                    end
                    alpha_d = alpha_d/prod(current_window);
                    alpha_f = alpha_f/prod(current_window);
                    % Compute variance
                    summation = 0;
                    for k=1:current_window(1)
                        for l=1:current_window(2)
                            if (alpha_f~=0)&&(alpha_d~=0)&&(sqrt((i-floor(window_size/2)+k-lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-lower_bound(2)-1)^2)~=0)
                                if l == 1
                                    summation = summation + ((right_window(k,2) - right_window(k,1))^2)/(alpha_f*alpha_d*sqrt((i-floor(window_size/2)+k-lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-lower_bound(2)-1)^2));
                                elseif l == current_window(2)
                                    summation = summation + ((right_window(k,l) - right_window(k,l-1))^2)/(alpha_f*alpha_d*sqrt((i-floor(window_size/2)+k-lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-lower_bound(2)-1)^2));
                                else
                                    summation = summation + (((right_window(k,l+1) - right_window(k,l-1))/2)^2)/(alpha_f*alpha_d*sqrt((i-floor(window_size/2)+k-lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-lower_bound(2)-1)^2));
                                end
                            end
                        end
                    end
                    % Check if the motion in the direction is possible
                    if (step==1)&&(summation==0)
                        break;
                    end
                    if step == 1
                        sigma = 1/summation;
                    elseif (1/summation) > sigma
                        motion_allowed(step) = 0;
                    else
                        var(step-1) = 1/summation;
                    end
                end
            end
            previous_window = start_window;
            previous_upper_bound = start_upper_bound;
            previous_lower_bound = start_lower_bound;
            if ~any(var~=Inf)
                if (step==1)&&(summation==0)
                    break;
                end
                continue;
            else
                min_var = find(var==min(var));
                for count=1:length(min_var)
                    % Update start_window
                    start_window = start_window + actions(min_var(count)+1,:);
                    % Update lower and upper bounds
                    if (min_var(count)==1)||(min_var(count)==2)
                        start_upper_bound = start_upper_bound + actions(min_var(count)+1,:);
                    elseif (min_var(count)==3)||(min_var(count)==4)
                        start_lower_bound = start_lower_bound - actions(min_var(count)+1,:);
                    end
                end
            end
        end
        if max(start_window)>window_size
            start_window = previous_window;
            start_upper_bound = previous_upper_bound;
            start_lower_bound = previous_lower_bound;
        end
        % Create the left window
        left_window = I_left_mod(i-start_lower_bound(1):i+start_upper_bound(1)-1,j-start_lower_bound(2):j+start_upper_bound(2)-1);
        % Create the right window
        right_window = I_right_mod(i-start_lower_bound(1):i+start_upper_bound(1)-1,j-start_lower_bound(2)+D_central:j+start_upper_bound(2)-1+D_central);
        % Compute alphas
        alpha_d = 0;
        alpha_f = 0;
        for k=1:start_window(1)
            for l=1:start_window(2)
                if l == 1
                    alpha_f = alpha_f + (right_window(k,2) - right_window(k,l))^2;
                elseif l == start_window(2)
                    alpha_f = alpha_f + (right_window(k,l) - right_window(k,l-1))^2;
                else
                    alpha_f = alpha_f + ((right_window(k,l+1) - right_window(k,l-1))/2)^2;
                end
                try
                    alpha_d = alpha_d + ((D(i-floor(window_size/2)+k-start_lower_bound(1)-1,j-floor(window_size/2)+l-start_lower_bound(2)-1) - D_central)^2)/sqrt((i-floor(window_size/2)+k-start_lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-start_lower_bound(2)-1)^2);
                catch
                    continue;
                end
            end
        end
        alpha_d = alpha_d/prod(start_window);
        alpha_f = alpha_f/prod(start_window);
        % Compute delta d
        numerator = 0;
        denominator = 0;
        for k=1:start_window(1)%-start_lower_bound(1):start_lower_bound(1)
            for l=1:start_window(2)%-start_lower_bound(2):start_lower_bound(2)
                if l == 1
                    numerator_num = (left_window(k,l) - right_window(k,l))*((right_window(k,l+1) - right_window(k,l))^2);
                    denominator_num = (right_window(k,l+1) - right_window(k,l))^2;
                elseif l == start_window(2)
                    numerator_num = (left_window(k,l) - right_window(k,l))*((right_window(k,l) - right_window(k,l-1))^2);
                    denominator_num = (right_window(k,l) - right_window(k,l-1))^2;
                else
                    numerator_num = (left_window(k,l) - right_window(k,l))*(((right_window(k,l+1) - right_window(k,l-1))/2)^2);
                    denominator_num = ((right_window(k,l+1) - right_window(k,l-1))/2)^2;
                end
                if (alpha_f~=0)&&(alpha_d~=0)&&(sqrt((i-floor(window_size/2)+k-start_lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-start_lower_bound(2)-1)^2)~=0)
                    numerator = numerator + numerator_num/(alpha_f*alpha_d*sqrt((i-floor(window_size/2)+k-start_lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-start_lower_bound(2)-1)^2));
                    denominator = denominator + denominator_num/(alpha_f*alpha_d*sqrt((i-floor(window_size/2)+k-start_lower_bound(1)-1)^2 + (j-floor(window_size/2)+l-start_lower_bound(2)-1)^2));
                end
            end
        end
        if (denominator ~= 0)&&(numerator~=0)
            D_improv(i-floor(window_size/2),j-floor(window_size/2)) = D_central + numerator/denominator;
        end
    end
end
D_improv = abs(D_improv);

% Change range to 0 to 255
D_improv = uint8(round(D_improv*255/max(max(D_improv))));

% Generate a histogram of the disparities
H = zeros(1,256);
for i=1:size(D_improv,1)
    for j=1:size(D_improv,2)
        H(D_improv(i,j)+1) = H(D_improv(i,j)+1) + 1;
    end
end

% Normalize the histogram to a value of 255
H = H*255/numel(D_improv);

% Generate the CDF of the histogram
new_H = zeros(1,256);
for i=1:256
    new_H(i) = round(sum(H(1:i)));
end

% Generate the new disparities using CDF
new_D = D_improv;
for i=1:size(new_D,1)
    for j=1:size(new_D,2)
        new_D(i,j) = new_H(D_improv(i,j)+1);
    end
end

% Display the image with colormap
imshow(new_D)
colormap(gca,jet)
colorbar

% Save the image
saveas(gca,['../output/SSD_adaptive_window/Disparity Map using Adaptive Window with SSD for window_size = ' num2str(window_size) '.jpg'])

end