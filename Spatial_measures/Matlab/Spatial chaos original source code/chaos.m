function [measure] = chaos(img,filter,varargin)
%% 0. Set environment variables
% addpath('Tools')
% pathName = fileparts(mfilename('fullpath'));
% addpath([pathName '\Tools'])

%% Parsing optional parameters
if (rem(length(varargin),2)==1)
    error('calcMoC: Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch varargin{i}
            case 'omega' % Parameter for bilateral filtering
                omega = varargin{i+1}; % Gaussian bilateral filter window
            case 'sigma_d' % Parameter for bilateral filtering
                sigma_d = varargin{i+1}; % Spatial-domain standard deviation
            case 'sigma_i' % Parameter for bilateral filtering
                sigma_i = varargin{i+1}; % Intensity-domain standard deviation
            case 'verbose' % Output option
                verbose = varargin{i+1};
            case 'quantVec' % Vector of quantile values with fixed stepsize 0.05
                quantVec = varargin{i+1};    
            case 'boundSkipFun' % Special function for skipping the boundary
                boundSkipFun = varargin{i+1}; % function handle
                if ~isa(varargin{i+1},'function_handle')
                   error('The boundSkip function needs to be a function handle');
                end 
            otherwise
                error(['calcMoC: Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

if ~exist('verbose','var'), verbose = 0; end

if ~exist('quantVec','var') % Vector of quantile values to take into account
    quantVec = 0.60:0.05:0.95;
end

if ~exist('omega','var') % Parameter for bilateral filtering
    omega = 200; % Gaussian bilateral filter window
    %omega = 10; % Gaussian bilateral filter window
end

if ~exist('sigma_d','var') % Parameter for bilateral filtering
    %sigma_d = 3; % Spatial-domain standard deviation
    sigma_d = 200; % Spatial-domain standard deviation
    
end

if ~exist('sigma_i','var') % Parameter for bilateral filtering
    %sigma_i = 0.3; % Intensity-domain standard deviation
    sigma_i = 0.5; % Intensity-domain standard deviation
end


%% 1. Apply edge detection on the image with prior bilateral filtering
img_save = img;

tic;
if filter == 1
    b = bilateral_colfilter(prepare_input(img), [omega omega], sigma_d, sigma_i); % bilateral filtering
elseif filter == 0
    b = prepare_input(img);
end

BF_time = toc;
tic;
[ax,ay] = gradient(b);
edge_crit = sqrt(ax.^2+ay.^2);
ED_time = toc;

edge_crit(isnan(img_save)) = NaN;
if exist('boundSkipFun','var')
    edge_crit = boundSkipFun(edge_crit);
end

%% 2. Now calculate the measure of chaos
ticID = tic;
[measRes,~,~,~] = calcMeasure(1,edge_crit,quantVec);
MoC_time = toc(ticID);
measure = measRes(1,1)-1;

if verbose
    disp(['Measure of chaos (MoC):' num2str(measure)]);
    disp('Needed');
    disp(['----------' num2str(BF_time+ED_time) 's for edge detection including bilateral filtering.']);
    disp(['----------' num2str(MoC_time) 's calculating the measure of chaos using the proposed method.']);
end


end

function c = prepare_input(img)
    img(isnan(img)) = 0; % Set all NaN's in the image as zeros
    c = img./(max(max(img))); % Normalize to [0,1]
    %c = abs(1-c);
end

function FA = bilateral_colfilter(A, size, std_c, std_s)

hs = (size-1)/2;
[x y] = meshgrid(-hs(2):hs(2),-hs(1):hs(1));
H_c = exp(-(x.*x + y.*y)/(2*std_c*std_c));
% perform filtering
FA = colfilt(A, size, 'sliding', @simfunc, std_s, H_c(:));
end

% adaptive similarity function
function V = simfunc(B, std, H_c);
% Find the row of B corresponding to the center pixel
center = floor((size(B,1)+1)/2);
% Make the similarity matrix
sim = exp(-(B - repmat(B(center,:), size(B,1), 1)).^2 / (2*std*std));
% Multiply the similarity matrix by the closeness matrix
sim = sim .* repmat(H_c, size(sim) ./ size(H_c));
ssum = sum(sim);
nonzerosum = ssum ~= 0;
sim(:, nonzerosum) = sim(:, nonzerosum) ./ repmat(ssum(nonzerosum), size(sim, 1), 1);
V = sum(sim .* B);
end