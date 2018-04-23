%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code recovers the model parameters for the defined number of
% Gaussian
% 
% Input:
%   data --> Data to be used to get the parameters
%      N --> Number of Gaussian
% 
% Output:
%         gmObj --> Gaussian Mixture Object. It contains the mean (N x D), 
%                  covariance (D x D x N), and mixture coefficients (N x 1)
%   isConverged --> Flag for checking ill-conditioned variance.
% 
% Submitted by: Ashwin Goyal (UID - 115526297)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gmObj, isConverged] = EM(data, N)

    % Ensure data is in the desired format
    if size(data,1) < size(data,2)
        data = data';
    end
    numSamples = size(data,1);
    % Initial mixture coefficients are equal
    mixtureCoeff = ones(N,1)/N;
    % Initial variances are equal
    for i = 1:N
        sigma(:,:,i) = diag(var(data));
    end
    % Initial mean is set as N random samples
    mu = datasample(data,N,1);
    
    % Start E-M
    oldLikelihood = -Inf;
    lastwarn('');
    while true
        %%% Estimation step
        % Compute all probabilities
        prob = zeros(numSamples,N);
        for i = 1:N
            prob(:,i) = mixtureCoeff(i)*gauss(data,mu(i,:),sigma(:,:,i));
            if any(prob(:,i) == Inf) || ~isempty(lastwarn)
                gmObj = [];
                isConverged = false;
                return
            end
        end
        % Compute posterior probability
        postProb = prob./sum(prob,2);
        
        %%% Verification
        % Compute likelihood
        likelihood = sum(log(sum(prob,2)));
        % Compute difference
        likelihoodDiff = likelihood - oldLikelihood;
        if (likelihoodDiff >= 0)&&(likelihoodDiff < abs(likelihood)*1e-6)
            isConverged = true;
            break;
        end
        oldLikelihood = likelihood;
        
        %%% Maximization step
        Nk = sum(postProb);
        for i = 1:N
            if Nk(i) == 0
                continue;
            end
            % Compute new mean
            mu(i,:) = sum(postProb(:,i).*data)/Nk(i);
            % Compute new variance
            sigma(:,:,i) = (sqrt(postProb(:,i)).*(data - mu(i,:)))'*(sqrt(postProb(:,i)).*(data - mu(i,:)))/Nk(i);
            % Compute new mixture coefficient
            mixtureCoeff(i) = Nk(i)/sum(Nk);
        end
    end
    gmObj = gmdistribution(mu,sigma,mixtureCoeff);

end

function N = gauss(x, mu, sigma)
% This function computes N(x|mu,sigma) for N-D data

    N = (1/(2*pi)^(size(x,2)/2))*(1/sqrt(det(sigma)))*exp(sum(-0.5*((x - mu)/sigma).*(x - mu),2));
    if any(isnan(N))
        N(isnan(N)) = Inf;
    end

end