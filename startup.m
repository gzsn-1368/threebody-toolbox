% Set environment variables.
format long;
toolboxRoot = fileparts(mfilename('fullpath'));
spiceKernelPath = fullfile(toolboxRoot, 'kernels');

requiredKernel = fullfile(spiceKernelPath, 'de440.bsp');
if ~isfile(requiredKernel)
    error('startup:MissingKernels', ...
        'Required SPICE kernels were not found in %s.', spiceKernelPath);
end

setenv('SPICE_KERNEL_PATH', spiceKernelPath);
