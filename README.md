# Threebody Toolbox Open-Source Release

This repository contains the open-source components of Threebody Toolbox.
The current release includes CRTBP dynamics, periodic-orbit initial-state
generation, differential correction, and the corresponding orbit data.

## Requirements

- MATLAB R2020b or newer is recommended.
- The examples in this release require MATLAB only.

## Installation

From MATLAB, run the installer in this repository:

```matlab
run('ToolboxInstall.m')
```

`ToolboxInstall.m` adds the repository root to the current MATLAB session and
does not modify the user's global MATLAB path. To make the path persistent,
use MATLAB's Set Path dialog or call `savepath` explicitly after reviewing the
path entry.

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

## Public API

| Function | Purpose |
| --- | --- |
| `dynamics.crtbp` | CRTBP state derivative, Jacobian, and variational equations |
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

## References

1. V. Szebehely, *Theory of Orbits: The Restricted Problem of Three Bodies*,
   Academic Press, 1967.
2. D. L. Richardson, “Analytic construction of periodic orbits about the
   collinear points,” *Celestial Mechanics*, 22, 241–253, 1980.
3. K. C. Howell, “Three-dimensional, periodic, halo orbits,” *Celestial
   Mechanics*, 32, 53–71, 1984.
4. W. S. Koon, M. W. Lo, J. E. Marsden, and S. D. Ross, *Dynamical Systems,
   the Three-Body Problem and Space Mission Design*, Springer, 2007.

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for contribution guidelines.
