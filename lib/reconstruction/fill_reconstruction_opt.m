function reconstruction_opt = fill_reconstruction_opt(reconstruction_opt)
if nargin<1
    reconstruction_opt = struct();
end

%% Method
reconstruction_opt.initial_learning_rate = ...
    default(reconstruction_opt,'initial_learning_rate',0.1);
reconstruction_opt.momentum = ...
    default(reconstruction_opt,'momentum',0.9);
reconstruction_opt.bold_driver_accelerator = ...
    default(reconstruction_opt,'bold_driver_accelerator',1.1);
reconstruction_opt.bold_driver_brake = ...
    default(reconstruction_opt,'bold_driver_brake',0.5);

%% Regularization
reconstruction_opt.is_regularized = ...
    default(reconstruction_opt,'is_regularized', ...
    isfield(reconstruction_opt,'regularizer'));
if reconstruction_opt.is_regularized
    reconstruction_opt.regularized = ...
        default(reconstruction_opt,'regularizer',0.1);
end

%% Soft thresholding
if isfield(reconstruction_opt,'target_min') && ...
        isfield(reconstruction_opt,'target_max')
    reconstruction_opt.is_thresholded = true;
    reconstruction_opt.soft_thresholding_factor = ...
        default(reconstruction_opt,'soft_thresholding_factor',0.1);
else
    reconstruction_opt.is_thresholded = false;
end

%% Verbosity
reconstruction_opt.is_verbose = default(reconstruction_opt,'is_verbose',true);
if reconstruction_opt.is_verbose && ...
        ~isfield(reconstruction_opt,'verbosity_period')
    reconstruction_opt.verbosity_period = 10;
end

%% Signal display
reconstruction_opt.is_signal_displayed = ...
    default(reconstruction_opt,'is_signal_displayed', ...
    isfield(reconstruction_opt,'signal_display_period'));
if reconstruction_opt.is_signal_displayed
    reconstruction_opt.signal_display_period = ...
        default(reconstruction_opt,'signal_display_period',10);
end

%% Alphanumeric ordering of field names
reconstruction_opt = orderfields(reconstruction_opt);
end
