# Threebody Toolbox Open-Source Release

This repository is the open-source release of the Threebody Toolbox.
It includes CRTBP orbit tools and a Moon-centred ephemeris model with its
force-model source, SPICE kernels, and gravity-field data.

## Requirements

- MATLAB R2020b or newer is recommended.
- A supported C compiler is required to build the ephemeris MEX file.
- Git LFS is required when cloning the repository because the SPICE kernels
  include large binary files.

## Installation

From MATLAB, run the installer in this repository:

```matlab
run('ToolboxInstall.m')
```

`ToolboxInstall.m` adds the repository root, the kernel directory, and the
CRTBP orbit-data directory to the MATLAB path. It runs `startup.m` to
configure `SPICE_KERNEL_PATH`, attempts to build the platform-specific
ephemeris MEX file if it is missing, and calls `savepath` so the paths remain
available in future MATLAB sessions.

The CRTBP functions can be used immediately. If automatic MEX compilation
fails, select a C compiler and retry the build manually:

```matlab
mex -setup C
run(fullfile('build', 'build_forcemodel.m'))
```

The build supports Apple Silicon and Intel macOS, 64-bit Linux, and 64-bit
Windows. On Windows, start MATLAB from an x64 Native Tools Command Prompt if
`cl.exe` and `lib.exe` are not already available on `PATH`.

## Quick start

```matlab
clear; clc; close all; format long;
run('ToolboxInstall.m');

odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[initialStates, periods] = orbits.crtbp_initstates( ...
    'halo', 'L2', 'northern', 15000, 0, ...
    'UsePlot', false, 'OdeOptions', odeOptions);

state0 = initialStates(1, :).';
[~, states] = ode113(@dynamics.crtbp, [0, periods(1)], state0, odeOptions);

plot3(states(:, 1), states(:, 2), states(:, 3), 'LineWidth', 1.0);
grid on; axis equal;
xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
title('Earth-Moon CRTBP halo orbit');
```

The same workflow is available in
[`examples/example_initial_state.m`](examples/example_initial_state.m) and
[`examples/example_differential_correction.m`](examples/example_differential_correction.m).

An ephemeris-model example is available in
[`examples/example_ephemeris.m`](examples/example_ephemeris.m).

## Public API

| Function | Purpose |
| --- | --- |
| `dynamics.crtbp` | CRTBP state derivative, Jacobian, and variational equations |
| `dynamics.ephem` | Moon-centred ephemeris dynamics and variational equations |
| `dynamics.ephemoptions` | Configure gravity, third-body, and radiation-pressure options |
| `orbits.crtbp_initstates` | Interpolate and correct periodic-orbit initial states |
| `orbits.crtbp_diffcorrection` | Correct a periodic-orbit state by single shooting |
| `constants` | Query physical and Earth--Moon CRTBP constants |

For detailed usage, run MATLAB's `help` command, for example
`help orbits.crtbp_initstates`, or read the help text at the beginning of the
corresponding function file.

The repository includes `resources/functionSignatures.json` for MATLAB
argument suggestions, making the public functions easier to discover and use.

## Units and orbit data

The CRTBP functions use the dimensionless Earth--Moon rotating frame. States
are ordered as:

```text
[x; y; z; dx/dt; dy/dt; dz/dt]
```

`orbits.crtbp_initstates` accepts amplitudes in kilometres and returns states
and periods in normalized CRTBP units. Conversion constants are available
through `constants`.

The `+orbits/crtbp_data` directory contains 41 files for 14 orbit families:

```text
axial, butterfly, dragonfly, dro, dro3d, dpo, halo,
longperiod, lpo, lyapunov, nrho, resonant, shortperiod, vertical
```

Each file contains `X`, `T`, `A`, `S`, `JacobiConstant`, and
`StabilityIndex`.

Some families have multiple branches for the same amplitude. When that
happens, the function returns multiple rows and records the matched values in
`info.targetamps` and `info.targetphases`.

`dynamics.ephem` uses SPICE ephemeris time in seconds past J2000 TDB. Position,
velocity, and acceleration are expressed in the Moon-centred J2000 inertial
frame in km, km/s, and km/s^2. Its kernels and gravity-field data are stored in
`kernels`; their paths are resolved relative to the repository.

### Supported force models

The ephemeris model always includes Earth and lunar gravity. The following
models can be configured with `dynamics.ephemoptions`:

- nonspherical Earth gravity using EGM2008, EGM96, or EIGEN-GL04C;
- nonspherical lunar gravity using LP100K or LP165P;
- gravitational perturbations from the Sun, Mercury, Venus, and the Mars,
  Jupiter, Saturn, Uranus, Neptune, and Pluto barycentres;
- relativistic acceleration correction;
- cannonball solar radiation pressure, including Earth and lunar eclipses;
- Earth-albedo acceleration.

By default, the model uses EGM2008 degree/order `[2, 0]`, LP165P degree/order
`[1, 0]`, solar gravity, and cannonball solar radiation pressure. Relativistic
and Earth-albedo corrections are disabled by default. A user-supplied control
acceleration can also be passed directly to `dynamics.ephem`.

## References

If you use this toolbox, cite at least one of References 1 and 2.

1. H. Hu, Z. Guo, Y. Liu, B. Jiang, G. Wen, and Y. Meng, “Relative motion
   configuration design and control of libration point orbits under
   high-fidelity ephemeris model,” *Aerospace Science and Technology*, 176,
   112697, 2026. https://doi.org/10.1016/j.ast.2026.112697
2. H. Hu, Z. Guo, Y. Liu, B. Jiang, G. Wen, and Y. Meng, “Fuel-optimal
   guidance between passively safe hold points for NRHO rendezvous,”
   *Aerospace Science and Technology*, 178, 113140, 2026.
   https://doi.org/10.1016/j.ast.2026.113140
3. V. Szebehely, *Theory of Orbits: The Restricted Problem of Three Bodies*,
   Academic Press, 1967.
4. D. L. Richardson, “Analytic construction of periodic orbits about the
   collinear points,” *Celestial Mechanics*, 22, 241–253, 1980.
5. K. C. Howell, “Three-dimensional, periodic, halo orbits,” *Celestial
   Mechanics*, 32, 53–71, 1984.
6. W. S. Koon, M. W. Lo, J. E. Marsden, and S. D. Ross, *Dynamical Systems,
   the Three-Body Problem and Space Mission Design*, Springer, 2007.
7. C. H. Acton, “Ancillary data services of NASA's Navigation and Ancillary
   Information Facility,” *Planetary and Space Science*, 44(1), 65–70, 1996.
8. C. H. Acton, N. J. Bachman, B. A. Semenov, and E. M. Wright, “A look
   toward the future in the handling of space science mission geometry,”
   *Planetary and Space Science*, 150, 9–12, 2018.

## License

The original Threebody Toolbox code is released under the MIT License. See
[`LICENSE`](LICENSE). CSPICE, kernels, and gravity-field data retain their
respective third-party terms and provenance; see
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for contribution guidelines.
