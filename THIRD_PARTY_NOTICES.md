# Third-party notices

The repository's MIT License applies only to the original Threebody Toolbox
code. The files listed below retain their original terms and provenance.

## CSPICE

`thirdparty/cspice` contains the source and public headers from CSPICE N0067,
published by NASA's Navigation and Ancillary Information Facility (NAIF). It
is included to support reproducible MEX builds and is not relicensed under the
Threebody Toolbox MIT License.

- Project: https://naif.jpl.nasa.gov/naif/toolkit.html
- Version: N0067
- Included version and package descriptions:
  `thirdparty/cspice/docs/version.txt` and
  `thirdparty/cspice/docs/dscriptn.txt`

Users should review NAIF's current usage rules and documentation before
redistributing CSPICE.

## SPICE kernels

The files below are distributed by NAIF and are used through CSPICE:

- `de440.bsp`
- `earth_1962_240827_2124_combined.bpc`
- `moon_pa_de440_200625.bpc`
- `moon_de440_220930.tf`
- `naif0012.tls`
- `pck00010.tpc`

The three `EM_*.fk`/`MoonCenteredInertial.fk` frame kernels define frames used
by Threebody Toolbox. SPICE kernel documentation and generic kernels are
available at https://naif.jpl.nasa.gov/pub/naif/generic_kernels/.

## Gravity-field data

The MATLAB data files `egm96.mat`, `egm2008.mat`, `eigengl04c.mat`,
`lp100k.mat`, and `lp165p.mat` contain coefficients converted from the named
Earth and lunar gravity-field models. The underlying scientific data retain
the attribution and terms of their originating agencies and model providers.
