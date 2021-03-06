function [yorig, kernel, deblur, opts] = ms_ngm_dirichlet_ubc_img(fn, opts)
% Do multi-scale blind deconvolution given input file name and options
% structure opts. Returns a double deblurred image along with estimated
% kernel. Following the kernel estimation, a non-blind deconvolution is run.
% The multi-scale framework is borrowed from Krishnan and Fergus CVPR 2011.
% Copyright (2015): Xu Zhou, Beihang University.
%

if (isempty(fn))
  if (isempty(opts.blur))
    fprintf('No image provided in fn or opts.blur!!!\n');
    return;
  else
    y = im2double(opts.blur);
  end;
else
  y = im2double(imread(fn));
end;


% save off for non-blind deconvolution
yorig = y;

%-------------------------------------------------------------------------
%                    parameters for kernel estimation
if ~isfield(opts.kernel_pars, 'lambda')
    opts.kernel_pars.lambda = 1e-6; % regularization weight on the negentropy term
end
if ~isfield(opts.kernel_pars, 'max_iter')
    opts.kernel_pars.max_iter = 20; % maximum number of iterations for alg.2
end
if ~isfield(opts.kernel_pars, 'back_alpha')
    opts.kernel_pars.back_alpha = 0.01; % back tracking line search
end
if ~isfield(opts.kernel_pars, 'back_beta') 
    opts.kernel_pars.back_beta = 0.5;% back tracking line search
end
if ~isfield(opts.kernel_pars, 'lower_bound') 
    opts.kernel_pars.lower_bound = 1; % decrease lb if the kernel is very large
end
if ~isfield(opts.kernel_pars, 'ng_min') 
    opts.kernel_pars.ng_min=1*10^(-5);% for stoppting criteria 
end
if ~isfield(opts.kernel_pars, 'cost_display') 
    opts.kernel_pars.cost_display=0;%  
end
if ~isfield(opts.kernel_pars, 'mode') 
    opts.kernel_pars.mode=0;% using the mode as kernel estimate?
end
if ~isfield(opts.kernel_pars, 'Laplacian_filter') 
    opts.kernel_pars.Laplacian_filter=[0 0 0;0 1 0;0 0 0];% use Gaussian prior on the kernel as default
end
%----------------------------------------------------------------------

%------------------------stop criteria---------------------------------
if ~isfield(opts,'k_tol')
    opts.k_tol=5e-4;
end
if ~isfield(opts,'xk_iter')
    opts.xk_iter = 20;
end


% gamma correct
y = y.^opts.gamma_correct;

% use a window to estimate kernel
if (~isempty(opts.kernel_est_win))
  w = opts.kernel_est_win;
  if (size(y, 3) == 3)
    y = rgb2gray(y(w(1):w(3), w(2):w(4), :));
  end;
else
  if (size(y, 3) == 3)
    y = rgb2gray(y);
  end;
end;

blur_size = opts.kernel_size;

% set kernel size for coarsest level - must be odd
% [ks1 ks2] = opt.kernel_size;
[max_ks, ind1] = max(opts.kernel_size);
if ind1 ==1 
    ind2 = 2;
else
    ind2 = 1;
end

minsize(ind1) = max(3, 2*floor(((max_ks - 1)/64)) + 1);
temp = floor(opts.kernel_size(ind2)/opts.kernel_size(ind1)*minsize(ind1));
if mod(temp,2)==0
    temp = temp +1;
end
minsize(ind2) = max(temp,3);

fprintf('Kernel size at coarsest level is [%d, %d]\n', minsize(1),minsize(2));

resize_step = sqrt(2);
% determine number of scales
num_scales = 1;
tmp = minsize(ind1);
while(tmp < max_ks)
  ksize(num_scales,ind1) = tmp;
  tmp2 = ceil(opts.kernel_size(ind2)/opts.kernel_size(ind1)*tmp);
  if mod(tmp2,2)==0
    tmp2 = tmp2 +1;
  end
  ksize(num_scales,ind2) = max(tmp2,3);
  
  num_scales = num_scales + 1;
  tmp = ceil(tmp * resize_step);
%     tmp = min(tmp+5,ceil(tmp * resize_step));
  if (mod(tmp, 2) == 0) 
    tmp = tmp + 1;
  end;
end;
ksize(num_scales,:) = opts.kernel_size;

ks = cell(num_scales);
alphas = ks;
ls = ks;

lambda_C = opts.kernel_pars.lambda_C; % increase it if the kernel is too sparse. 
% blind deconvolution - multiscale processing
for s = 1:num_scales
  if (s == 1)
    % at coarsest level, initialize kernel
    if max_ks>50
        Gsigma = 1;
    else
        Gsigma = 0.5;
    end
    ks{s} = init_kernel(ksize(1,:),Gsigma);
    alphas{s} = ks{s}+opts.kernel_pars.lower_bound;
    k1 = ksize(1,1);
    k2 = ksize(1,2); 
  else
    % upsample kernel from previous level to next finer level
    k1 = ksize(s,1);
    k2 = ksize(s,2); 
    
    % resize kernel from previous level
    tmp = ks{s-1};
    tmp(tmp<0) = 0;
    tmp = tmp/sum(tmp(:));
    ks{s} = imresize(tmp, [k1 k2], 'bilinear');
    alphas{s} = imresize(alphas{s-1}, [k1 k2], 'bilinear');
    % bilinear interpolantion not guaranteed to sum to 1 - so renormalize
    ks{s}(ks{s} < 0) = 0;
    sumk = sum(ks{s}(:));
    ks{s} = ks{s}./sumk;
  end;
  
  % image size at this level
  r = floor(size(y, 1) * k1 / blur_size(1));
  c = floor(size(y, 2) * k2 / blur_size(2));
  
  if (s == num_scales)
    r = size(y, 1);
    c = size(y, 2);
  end;
  
  fprintf('Processing scale %d/%d; kernel size %dx%d; image size %dx%d\n', ...
            s, num_scales, k1, k2, r, c);
  
  % resize y according to the ratio of filter sizes
  ys = imresize(y, [r c], 'bilinear');

  if (s == 1)
    ls{s} = ys;
  else
      % upscale the estimated derivative image from previous level
      ls{s} = imresize(ls{s - 1}, [r c], 'bilinear');

  end
  if (s == num_scales)
      opts.kernel_pars.lambda_C=lambda_C;
  else
      % update the weight lambda_h according to the size of kernel or image
      opts.kernel_pars.lambda_C=lambda_C*ksize(s,1)*ksize(s,2)/(ksize(num_scales,1)*ksize(num_scales,2));
%       opts.kernel_pars.lambda_C=lambda_C*r*c/numel(y);
  end
  
     % center the kernel
  [ls{s}, ks{s}, shift_kernel] = center_kernel_img_space(ls{s},ks{s});
  alphas{s}=max(conv2(alphas{s},shift_kernel,'same'),opts.kernel_pars.lower_bound);
  
  % call kernel estimation for this scale
   opts.kernel_pars.alpha0=alphas{s};
  [ls{s} ks{s} alphas{s}] = ss_ngm_dirichlet_ubc_img(ys, ls{s}, ks{s},alphas{s},opts);
  if (s == num_scales)
    kernel = alphas{s}-opts.kernel_pars.lower_bound;    
    kernel = kernel / sum(kernel(:));
  end
end
if nargout>=3
    deblur = yorig;
    if opts.use_ycbcr
      if (size(yorig, 3) == 3)
        ycbcr = rgb2ycbcr(yorig);
      else
        ycbcr = yorig;
      end
    end
    if opts.use_ycbcr
       [x_fov] = firls_deb_ubc(ycbcr(:,:,1),kernel,opts.firls);    
        deblur(:, :, 1) = x_fov;
        if (size(ycbcr, 3) == 3)
          deblur(:, :, 2:3) = ycbcr(:, :, 2:3);
          deblur = ycbcr2rgb(deblur);
        end
    else
      for j = 1:size(yorig,3)
        [x_fov] = firls_deb_ubc(yorig(:,:,j),kernel,opts.firls);
        deblur(:, :, j) = x_fov;
      end
    end
end
opts.ls=ls{s};
end

function [k] = init_kernel(minsize,Gsigma)
%    k = zeros(minsize, minsize);
%    k((minsize - 1)/2, (minsize - 1)/2:(minsize - 1)/2+1) = 1/2;
%    k0 = ones(3)/9;
%    k = zeros(minsize);
%    k((minsize(1) - 1)/2:(minsize(1) - 1)/2+2, (minsize(2) - 1)/2:(minsize(2) - 1)/2+2) = k0;
  k=fspecial('gaussian',minsize,Gsigma);

end
