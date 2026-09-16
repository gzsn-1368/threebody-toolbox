# Contributing

Thank you for helping improve the Threebody Toolbox open-source release.

Keep contributions focused on reproducible CRTBP orbit generation,
differential correction, orbit data, or documentation. Do not add local
absolute paths, proprietary data, compiled binaries, SPICE kernels, or
unreviewed third-party code.

For pull requests:

- Explain the numerical or documentation change in the summary.
- Include a minimal MATLAB reproduction for behavior changes.
- State the MATLAB release and additional toolboxes used for verification.
- Report tolerances, residuals, and any changes to units or coordinate
  conventions.
- Update the README or function help when the public interface changes.

Use English help text, portable paths constructed with `fullfile`, descriptive
variable names, and deterministic examples. Preserve the normalized
Earth--Moon rotating-frame convention unless a change is explicitly
documented.
