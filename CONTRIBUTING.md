# Contributing

Thank you for helping improve the Threebody Toolbox open-source release.

Keep contributions focused on reproducible astrodynamics models, orbit
generation, correction methods, data, or documentation. Do not add local
absolute paths, proprietary data, compiled binaries, or unreviewed
third-party code. Kernel or third-party updates must include their source,
version, and redistribution notice.

For pull requests:

- Explain the numerical or documentation change in the summary.
- Include a minimal MATLAB reproduction for behavior changes.
- State the MATLAB release and additional toolboxes used for verification.
- Report tolerances, residuals, and any changes to units or coordinate
  conventions.
- Update the README or function help when the public interface changes.

Use English help text, portable paths constructed with `fullfile`, descriptive
variable names, and deterministic examples. State the coordinate frame and
units for every public dynamics interface.
