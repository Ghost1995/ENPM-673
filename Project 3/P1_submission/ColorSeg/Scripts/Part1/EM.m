%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code recovers the model parameters for the defined number of
% gaussians
% 
% Input:
%           N --> Number of Gaussians
%        data --> Data to be used to get the parameters
%   plot_path --> The relative path to save the figure
% 
% Input:
%      mu --> Mean of the Gaussians (N x D)
%   sigma --> Covariance of the Gaussians (N x D x D)
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [mu, sigma] = EM(N, data, plot_path)

N = 5;
plot_path = '';
mu = rand(N,3)*1;
sigma = rand(N,3,3)*1;
data1(:,1) = normrnd(mu(1,1),sigma(1,1,1),[10 1]);
data1(:,2) = normrnd(mu(1,2),sigma(1,2,2),[10 1]);
data1(:,3) = normrnd(mu(1,3),sigma(1,3,3),[10 1]);
data2(:,1) = normrnd(mu(2,1),sigma(2,1,1),[10 1]);
data2(:,2) = normrnd(mu(2,2),sigma(2,2,2),[10 1]);
data2(:,3) = normrnd(mu(2,3),sigma(2,3,3),[10 1]);
data3(:,1) = normrnd(mu(3,1),sigma(3,1,1),[10 1]);
data3(:,2) = normrnd(mu(3,2),sigma(3,2,2),[10 1]);
data3(:,3) = normrnd(mu(3,3),sigma(3,3,3),[10 1]);
% data1 = normrnd(reshape(repmat(mu(1,:),2,1),[1,2,2]),sigma(1,:,:),[10 2]);
% data2 = normrnd(reshape(repmat(mu(2,:),2,1),[1,2,2]),sigma(2,:,:),[10 2]);
% data3 = normrnd(reshape(repmat(mu(3,:),2,1),[1,2,2]),sigma(3,:,:),[10 2]);
data = [data1;data2;data3];
init_mu = mu;
init_sigma = sigma;
% plot(data1,data1,'r*')
% hold on
% viscircles([mu(1) mu(1)],sigma(1),'Color','r')
% plot(data2,data2,'g*')
% viscircles([mu(2) mu(2)],sigma(2),'Color','g')
% plot(data3,data3,'b*')
% viscircles([mu(3) mu(3)],sigma(3),'Color','b')
% hold off
    
    if size(data,1) < size(data,2)
        data = data';
    end
    mu = rand(N,size(data,2)).*max(data);
    sigma = sqrt(var(data));%zeros(N,size(data,2),size(data,2));
%     for i = 1:N
%         for j = 1:size(data,1)
%             sigma(i,:,:) = sigma(i,:,:) + reshape((data(j,:) - mu(i,:))'*(data(j,:) - mu(i,:)),[1,size(data,2),size(data,2)]);
%         end
%     end
%     sigma = sigma/size(data,1);
    % mu = mean(data)*ones(N,size(data,2));
    % sigma = cov(data)*ones(N,size(data,2),size(data,2));
    % mu(1) = mu(2) - sigma(1);
    % mu(3) = mu(2) + sigma(3);
    mixtureCoeff = ones(N,1)/N;
    % pk = pk/sum(pk);
    likelihood = 0;
    for i = 1:size(data,1)
        indLikelihood = 0;
        for j = 1:N
            indLikelihood = indLikelihood + mixtureCoeff(j)*gauss(data(i,:),mu(j,:),sigma(j,:,:));
        end
        likelihood = likelihood + log(indLikelihood);
    end
    old_likelihood = -Inf;
    
    while (likelihood > old_likelihood)||(abs(likelihood) == Inf)
        old_mu = mu;
        old_sigma = sigma;
        old_likelihood = likelihood;
    %     plot(data,data,'*')
    %     hold on
    %     viscircles([mu' mu'],sigma);
    %     hold off
        prob = zeros(size(data,1),N);
        for i = 1:size(data,1)
            indProb = zeros(N,1);
            for j = 1:N
                indProb(j) = mixtureCoeff(j)*gauss(data(i,:),mu(j,:),sigma(j,:,:));
            end
            if (sum(indProb) ~= 0)%&&(abs(indProb) ~= Inf)
                prob(i,:) = (mixtureCoeff.*indProb)'/sum(indProb);
            end
        end
        num = sum(prob);
        num(num == 0) = (size(data,1) - sum(num))/sum(num == 0);
        for i = 1:N
            mu(i,:) = sum(repmat(prob(:,i),1,size(data,2)).*data)/num(i);
            temp_sigma = zeros(size(data,2), size(data,2));
            for j = 1:size(data,1)
                temp_sigma = temp_sigma + prob(j,i)*(data(j,:) - mu(i,:))'*(data(j,:) - mu(i,:));
            end
            sigma(i,:,:) = reshape(temp_sigma/num(i),[1,size(temp_sigma)]);
            mixtureCoeff(i) = num(i)/size(data,1);
        end
        likelihood = 0;
        for i = 1:size(data,1)
            indLikelihood = 0;
            for j = 1:N
                indLikelihood = indLikelihood + mixtureCoeff(j)*gauss(data(i,:),mu(j,:),sigma(j,:,:));
            end
            likelihood = likelihood + log(indLikelihood);
        end
    end
    mu = old_mu;
    sigma = old_sigma;
    if (~strcmp(plot_path,''))||(~isempty(plot_path))
        figure
        hold on
        for i = 1:N
            for j = 1:size(data,2)
                data(:,j) = normrnd(mu(i,j),sigma(i,j,j),[10000 1]);
            end
            if size(data,2) == 1
                plot(data)
            elseif size(data,2) == 2
                plot(data(:,1),data(:,2))
            elseif size(data,2) == 3
                plot3(data(:,1),data(:,2),data(:,3))
            end
        end
        hold off
        title('Computed Gaussians')
        saveas(gcf,[plot_path 'EM' num2str(size(data,2)) 'D' num2str(N) 'N.jpg'])
    end
% end


function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma)

    sigma = reshape(sigma,[size(sigma,2) size(sigma,3)]);
    N = (1/(2*pi)^(length(x)/2))*(1/sqrt(det(sigma)))*exp(-0.5*((x - mu)/sigma)*(x - mu)');
    if isnan(N)
        N = 0;
    end
    
end