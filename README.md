# MC Dataset

MC Dataset is an extensive collection of datasets computed with
[PyXOpto](https://github.com/xopto/pyxopto) Monte Carlo light propagation
models. The readily available datasets are computed for a vast variety of
sources, detectors and sample optical properties and include information of
reflectance, transmittance, sampling volume and fluence/energy deposition data.
The datasets can be easily customized through the 
[dataset](https://xopto.github.io/pyxopto/docs/html/source/dataset/index.html)
module of the [PyXOpto](https://github.com/xopto/pyxopto) project.

## Links to the data files of the latest release

Datasets are split into multiple files to reduce the download size and to allow
partial downloads. The following table has the links to the dataset files of
the latest release. A more detailed description of the individual datasets is
given in the following section.

| Dataset                                                       | Release file                                                                                                                                                                          | Description                                                                                  |
|---------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
|  [MCML comparison datasets](#mcml-comparison-datasets)        |                                                                                                                                                                                       |                                                                                              |
|                                                               | [mcml_comparison.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_comparison.zip)                                                                                | single layer 100 mm and 1 mm thick, 2-layers 0.1 mm and 1 mm thick                           |
|  [Layered media datasets (MCML)](#layered-media-datasets)     |                                                                                                                                                                                       |                                                                                              |
|                                                               | [mcml_1-layer-semiinfinite_line.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_line.zip)                                                  | Semi-infinite medium, line source                                                            |
|                                                               | [mcml_1-layer-semiinfinite_collimated-200um.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_collimated-200um.zip)                          | Semi-infinite medium, collimated source with 200 µm diameter                                 |
|                                                               | [mcml_1-layer-semiinfinite_gaussian-fwhm-100um.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_gaussian-fwhm-100um.zip)                    | Semi-infinite medium, Gaussian source with 100 µm FWH                                        |
|                                                               | [mcml_1-layer-semiinfinite_fiber-200um-0.22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_fiber-200um-0.22na.zip)                      | Semi-infinite medium, fiber source with 200 µm 0.22 NA core                                  |
|                                                               | [mcml_1-layer-semiinfinite_six-linear-200um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_six-linear-200um-0_22na.zip)            | Semi-infinite medium, six-linear array probe, 200 µm 0.22 NA fiber core                      |
|                                                               | [mcml_1-layer-semiinfinite_six-around-one-200um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_six-around-one-200um-0_22na.zip)    | Semi-infinite medium, six-arround-one probe, 200 µm 0.22 NA fiber core                       |
|                                                               | [mcml_1-layer-semiinfinite_six-around-one-400um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_six-around-one-400um-0_22na.zip)    | Semi-infinite medium, six-arround-one probe, 400 µm 0.22 NA fiber core                       |
|                                                               | [mcml_1-layer-semiinfinite_single-fiber-100um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_single-fiber-100um-0_22na.zip)        | Semi-infinite medium, single-fiber probe, 100 µm 0.22 NA fiber core                          |
|                                                               | [mcml_1-layer-semiinfinite_single-fiber-200um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_single-fiber-200um-0_22na.zip)        | Semi-infinite medium, single-fiber probe, 200 µm 0.22 NA fiber core                          |
|                                                               | [mcml_1-layer-semiinfinite_single-fiber-400um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_single-fiber-400um-0_22na.zip)        | Semi-infinite medium, single-fiber probe, 400 µm 0.22 NA fiber core                          |
|                                                               | [mcml_1-layer-semiinfinite_single-fiber-800um-0_22na.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcml_1-layer-semiinfinite_single-fiber-800um-0_22na.zip)        | Semi-infinite medium, single-fiber probe, 800 µm 0.22 NA fiber core                          |
| [Voxelized media datasets (MCVOX)](#voxelized-media-datasets) |                                                                                                                                                                                       |                                                                                              |
|                                                               | [mcvox.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcvox.zip)                                                                                                    | Energy deposition for a voxelized 2-layer skin with an embedded blood vessel at 26 depths]() |                                                                                              |
|                                                               | [mcvox-2-layer-skin-200um-vessel-500um-depth-deposition.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/mcvox-2-layer-skin-200um-vessel-500um-depth-deposition.zip)  | Energy deposition for a voxelized 2-layer skin with an embedded blood vessel at 0.5 mm depth |
| [Sampling volume datasets (SV)](#sampling-volume-datasets)    |                                                                                                                                                                                       |                                                                                              |
|                                                               | [sv.zip](https://github.com/xopto/mcdataset/releases/download/v0.0.1/sv.zip)                                                                                                          | Sampling volume dataset                                                                      |

# Datasets

The following sections provide information on the  four groups of datasets:
 * [MCML comparison datasets](#mcml-comparison-datasets)
 * [Layered media datasets (MCML)](#layered-media-datasets)
 * [Voxelized media datasets (MCVOX)](#voxelized-media-datasets)
 * [Sampling volume datasets (SV)](#sampling-volume-datasets)

At the end of each section there is a breakdown of the directory
tree and naming conventions that are used to organize the dataset files. 

## MCML comparison datasets

This dataset is intended for comparison with the
[MCML](https://omlc.org/software/mc/mcml/index.html) package. Datasets are
computed for a 1-layer (1&nbsp;mm and 100&nbsp;mm thick) and 2-layer sample
(0.1&nbsp;mm and 1.0&nbsp;mm thick layers). The optical properties are varied
in the top sample layer:
* Absorption coefficient is sampled from [0.0, 2.5, 5.0]&nbsp;cm <sup>-1</sup>.
* Reduced scattering coefficient is sampled from [5.0, 20.0, 35.0]&nbsp;cm <sup>-1</sup>.
* Scattering phase function anisotropy is sampled from [0.1, 0.5, 0.9].
* The refractive index of the sample is set to:
    * 1-layer sample: 1.337,
    * 2-layer sample: 1.462 for the top layer and 1.337 for the bottom layer.
* The refractive index of the surrounding medium is set to 1.0.
* The remaining optical properties of the 2<sup>nd</sup> layer are fixed
to *µ<sub>a</sub>*=0.5&nbsp;cm<sup>-1</sup>,
*µ<sub>s</sub>'*=20.0&nbsp;cm<sup>-1</sup> and *g*=0.8.

The datasets are computed with a normally incident infinitely thin source. The
reflectance and transmittance are collected through radial detectors with 500
concentric accumulators from 0 to 5&nbsp;mm.

![Radial detector](/docs/source/static/mcdetector_radial_start_0.svg)<br>
*Illustration of the radial detector that is used to collect the reflectance and transmittance as `Radial(Axis(start=0.0, stop=0.005, n=500))`. The central and outermost accumulator are highlighted with a gray fill*

The energy deposition is collected
in a 2D accumulator array along the radial *r* and *z* axis. The range of both
axis is from 0 to 5&nbsp;mm and the number of accumulators along each axis
is set 500. All the simulations are prepared with 100 million photon packets
and simulation termination radius is set to 1000&nbsp;mm.

![Deposition detector](/docs/source/static/mcfluence_fluence_rz_3d.svg)<br>
*Illustration of the radially symmetric energy deposition detector `FluenceRz`.*

### Dataset files

Datasets are available as compressed numpy data files that are organized as
follows:

```
data/mcml_comparison/<sample>/line/radial/hg/g-<g>/mua-<mua>-musr-<musr>-invcm.npz
```

The values of placeholders <> are as follows:

* `<sample>` can take one of the following values:

    * `1-layer-1mm` 

        A single layer 1&nbsp;mm thick medium.

    * `1-layer-100mm`

      A single layer 100&nbsp;mm thick medium.

    * `2-layer-100um-1mm`

      A two-layer medium with 0.1&nbsp;mm thick top layer and 1&nbsp;mm thick bottom layer

* `<g>` is the anisotropy formatted with two decimal digits and `_`
  as the decimal separator, e.g `0_10` for *g*=0.1, 
  `0_50` for *g*=0.5 or `0_90` for *g*=0.9.

* `<mua>` is the absorption coefficient in units of cm&nbsp;<sup>-1</sup> with
  formatted with two decimal digits and `_` as the decimal separator, e.g `2_50`
  for *µ<sub>a</sub>*=2.5&nbsp;cm<sup>-1</sup>.

* `<musr>` is the reduced scattering coefficient in units of cm&nbsp;<sup>-1</sup> with
  formatted with two decimal digits and `_` as the decimal separator, e.g `25_00`
  for *µ<sub>s</sub>'*=25.0&nbsp;cm<sup>-1</sup>.

![HG phase function](/docs/source/static/mcml_comparison_hg.svg)<br>
*Henyey-Greenstein scattering phase functions used in the MCML comparison dataset.*

## Layered media datasets (MCML)

The datasets are produced with 100 million photon packets, except for the
Spatial Frequency Domain Imaging (SFDI) dataset and all the datasets that use an optical
fiber probe with a linear layout of 6 fibers. These
datasets are run with 1000 million photon packets. The simulation termination
radius is set to 25&nbsp;mm, except for the spatial frequency domain imaging (SFDI)
dataset that uses a 150&nbsp;mm simulation termination radius. The optical
properties of the sample are varied according to the values listed in the following
two tables.

*Absorption and reduced scattering coefficients*
| Parameter        | From (cm<sup>-1</sup>) | To (cm<sup>-1</sup>) | Points |
|------------------|------------------------|----------------------|--------|
| *µ<sub>a</sub>*  | 0.0                    | 5.0                  | 21     |
| *µ<sub>s</sub>'* | 5.0                    | 35.0                 | 21     |

*Scattering phase functions*

| Scattering phase function        | Parameter              | Values                          |
|----------------------------------|------------------------|---------------------------------|
| HG (Henyey-Greenstein)           | *g*                    | 0.1, 0.3, 0.5, 0.7, 0.9         |
| MHG (Modified Henyey-Greenstein) | *g*                    | 0.1, 0.3, 0.5, 0.7, 0.9         |
|                                  | *β*                    | 0.0, 0.2, 0.4, 0.6, 0.8, 1.0    |
| GK (Gegenbauer Kernel)           | *g*                    | 0.1, 0.3, 0.5, 0.7, 0.9         |
|                                  | *α*                    | -0.5, 0.0, 0.5, 1.0, 1.5        |
| MIE-polystyrene                  | *λ*                    | 500&nbsp;nm                     |
|                                  | *diameter*             | 0.25, 0.5 1.0, 2.0, 4.0&nbsp;µm |
|                                  | *n<sub>particle</sub>* | 1.603                           |
|                                  | *n<sub>medium</sub>*   | 1.337                           |
| MIE-fused silica                 | *λ*                    | 500&nbsp;nm                     |
|                                  | *diameter*             | 0.25, 0.5 1.0, 2.0, 4.0&nbsp;µm |
|                                  | *n<sub>particle</sub>* | 1.462                           |
|                                  | *n<sub>medium</sub>*   | 1.337                           |

The refractive index of the sample is set to 1.337.

![HG phase function](/docs/source/static/mcml_hg.svg)<br>
*Henyey-Greenstein (HG) scattering phase functions as defined in the above table.*

![MHG phase function](/docs/source/static/mcml_mhg.svg)<br>
*Examples of Modified Henyey-Greenstein (MHG) scattering phase functions from the above table.*

![GK phase function](/docs/source/static/mcml_gk.svg)<br>
*Examples of Gegenbauer Kernel (GK) scattering phase functions from the above table.*

![Mie fused silica phase function](/docs/source/static/mcml_mie_fused_silica.svg)<br>
*Mie scattering phase functions of water-suspended fused silica particles specified in the above table.*

![Mie fused silica phase function](/docs/source/static/mcml_mie_polystyrene.svg)<br>
*Mie scattering phase functions of water-suspended polystyrene particles specified in the above table.*

Datasets are available for the following basic sources that use a laterally
uniform boundary between the sample and the surrounding medium.

*Basic sources with a uniform sample-source interface and related reflectance detectors*
| Source       | Parameter   | Value       | Reflectance detector        |
|--------------|-------------|-------------|-----------------------------|
| Line         |             |             | Radial(Axis(0, 0.005, 500)) |
| UniformBeam  | *diameter*  | 200&nbsp;μm | Radial(Axis(0, 0.005, 500)) |
| GaussianBeam | *FWHM*      | 100&nbsp;μm | Radial(Axis(0, 0.005, 500)) |
| UniformFiber | *dcore*     | 200&nbsp;μm | Radial(Axis(0, 0.005, 500)) |
|              | *dcladding* | 220&nbsp;μm |                             |
|              | *ncore*     | 1.462       |                             |
|              | *NA*        | 0.22        |                             |

The refractive index of the surrounding medium is set to 1.0 except when
using the UniformFiber source,
when the refractive index of the surrounding medium follows the
refractive index of the fiber core 1.462.

The reflectance of basic sources is collected with a radial detector
with range from 0 to 5&nbsp;mm and 500 concentric accumulators, each 5&nbsp;μm
wide. The acceptance angle is unlimited, except for the 
UniformFiber source for which it is limited by the NA of the fiber core. The acceptance angle within the fiber core.

Datasets are also prepared for optical fiber probe sources that use a surface
layout to more accurately describe the interface between the optical fiber
probe tip and the sample. All the probe sources launch the photon packets with
the UniformFiber source.

*Optical fiber probe sources with a detailed sample-source interface and*
*related reflectance detectors*
| Probe         | Parameter       | Value        | Description               | Reflectance detector |
|---------------|-----------------|--------------|---------------------------|----------------------|
| SixAroundOne  | *dcore*         | 200&nbsp;μm  | six-around-one layout     | SixAroundOne         |
|               | *dcladding*     | 220&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| SixAroundOne  | *dcore*         | 400&nbsp;μm  | six-around-one layout     | SixAroundOne         |
|               | *dcladding*     | 420&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| LinearArray   | *dcore*         | 200&nbsp;μm  | linear layout of 6 fibers | LinearArray          |
|               | *dcladding*     | 220&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *n*             | 6            |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| LinearArray   | *dcore*         | 100&nbsp;μm  | single fiber layout       | LinearArray          |
|               | *dcladding*     | 120&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *n*             | 1            |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| LinearArray   | *dcore*         | 200&nbsp;μm  | single fiber layout       | LinearArray          |
|               | *dcladding*     | 220&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *n*             | 1            |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| LinearArray   | *dcore*         | 400&nbsp;μm  | single fiber layout       | LinearArray          |
|               | *dcladding*     | 420&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *n*             | 1            |                           |                      |
|               | *diameter*      | 6.0&nbsp;mm  |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |
| LinearArray   | *dcore*         | 800&nbsp;μm  | single fiber layout       | LinearArray          |
|               | *dcladding*     | 820&nbsp;μm  |                           |                      |
|               | *ncore*         | 1.462        |                           |                      |
|               | *NA*            | 0.22         |                           |                      |
|               | *n*             | 1            |                           |                      |
|               | *diameter*      | 6.0&nbsp; mm |                           |                      |
|               | *reflectivity*  | 0.6          |                           |                      |

![Six-around-one layout](/docs/source/static/mcsurface_six-around-one-200-400um.svg)<br>
*Illustration of the surface layouts for the two six-around-one optical fiber probes (`SixAroundOne`) with 200 and 400&nbsp;μm fiber cores.*

![Linear array of 6 fibers](/docs/source/static/mcsurface_six-linear-200um.svg)<br>
*Illustration of the surface layout for the linear array probe (`LinearArray`) with 6 optical fibers.*

![Linear array of 6 fibers](/docs/source/static/mcsurface_single-fiber-100-200-400-800um.svg)<br>
*Illustration of the surface layouts for the single fiber probes (`LinearArray`) with 100, 200, 400 and 800&nbsp;μm fiber cores.*

The reflectance of optical fiber probe sources is collected only through
the individual optical fibers of the probe.

### SFDI source-detector arrangements

The SFDI datasets are computed for two source-detector configurations and
include raw reflectance and the corresponding frequency-domain reflectance,
which is computed for spatial frequencies from 0.00 to 0.80&nbsp;mm<sup>-1</sup>
in 0.01&nbsp;mm<sup>-1</sup> steps:

* A normally incident Line source and a radial 
  `Radial(Axis(0.0, stop=0.15, n=4000, logscale=True), cosmin=0.98481)`
  reflectance detector that uses 4000 logarithmically spaced concentric
  accumulators from 0 to 150&nbsp;mm. The acceptance angle is limited to 10°.
  Hankel transform is used to compute the spatial frequency-domain reflectance.
  Note that this transform produces real values. 

![Sampling volume](/docs/source/static/mcdetector_radial_logscale_start_0.svg)<br>
*A Radial detector with logarithmically spaced accumulators. The outermost accumulator is highlighted with a gray fill*

* A normally incident Line source and a tilted linear detector with 20°
  incidence (along the *x* axis). The accumulators of the detector extend to
  infinity along the positive and negative *y* axis and follow a logarithmic
  spacing along the positive and negative direction of the *x* axis
  `SymmetricX(SymmetricAxis(center=0.0, range=0.15, n_half=4000, logscale=True), cosmin=0.98480)`.
  The described detector uses 8000 (4000 in each direction along the
  :math:`x` axis) logarithmically spaced accumulators *x*=-150
  to *x*=150;&nbsp;mm. The acceptance angle of the detector is limited to 10°
  around the tilted detector axis. Fourier transform is used to compute
  the spatial frequency-domain reflectance. Note that this transform
  produces complex values with amplitude and phase information.

![Sampling volume](/docs/source/static/mcdetector_symmetricx_logscale.svg)<br>
*A SymmetricX detector with logarithmically spaced accumulators around the central axis. The first and last accumulators are highlighted with a gray fill.*

Note that the SFDI datasets are run with 1000 million photon packets and that
the simulation termination radius is set to 150&nbsp;mm.

### Dataset files

Datasets are available as compressed numpy data files that are organized as
follows:

```
data/mcml/<sample>/<source>/<detector>/<pf>/<pf_param_1>/<pf_param_2>/mua-<mua>-musr-<musr>-invcm.npz
```

The values of placeholders <> are as follows:

* `<sample>` can take the following values:

    * `1-layer-semiinfinite`

        A single sample layer of infinite thickness.

* `<source>` is the type of the photon packet source / probe used in the datasets:

    * `line`

       Infinitely narrow line source (`Line`).

    * `collimated-200um`

        A 200&nbsp; µm beam diameter (UniformBeam).

    * `gaussian-fwhm-100um`

        A Gaussian beam with 100&nbsp;µm FWHM (`GaussianBeam`).

    * `fiber-200um-0_22na`

        A 200&nbsp;µm, 0.22 NA fiber source (`UniformFiber`).

    * `six-around-one-200um-0_22na` 

        A six-around-one layout optical fibers
        (200&nbsp;µm core, 220&nbsp;µm cladding, 0.22 NA) (`SixAroundOne`)
        with the central optical fiber used as the source (`UniformFiber`).
        The fibers are tightly packed.

    * `six-around-one-400um-0_22na`

        For a six-around-one layout optical fibers
        (400&nbsp;µm core, 420&nbsp;µm cladding, 0.22 NA) (`SixAroundOne`)
        with the central optical fiber used as the source (`UniformFiber`).
        The fibers are tightly packed.

    * `six-linear-array-200um-0_22na`

        For a linear layout of 6 optical fibers
        (200&nbsp;µm core, 220&nbsp;µm cladding, 0.22 NA) (`LinearArray`)
        with the leftmost optical fiber used as the source (`UniformFiber`).
        The fibers are tightly packed.

    * `single-fiber-100um-0_22na`

        A single fiber layout
        (100&nbsp;µm core, 120&nbsp;µm cladding, 0.22 NA) (`LinearArray`).

    * `single-fiber-200um-0_22na`

        A single fiber layout
        (200&nbsp;µm core, 220&nbsp;µm cladding, 0.22 NA) (`LinearArray`).

    * `single-fiber-400um-0_22na`

        A single fiber layout
        (400&nbsp;µm core, 420&nbsp;µm cladding, 0.22 NA) (`LinearArray`).

    * `single-fiber-400um-0_22na`

        A single fiber layout
        (800&nbsp;µm core, 820&nbsp;µm cladding, 0.22 NA) (`LinearArray`).

* `<detector>` is the type of detector used by the datasets:
    
    * `radial`

        For simple sources with laterally uniform source-sample boundary,

    * `probe`

        For optical fiber probes with surface layout.

* `<pf>` is the type of scattering phase function used in the datasets:
    
    * `hg` for HG.

    * `mhg` for MHG.

    * `gk` for GK.

    * `mie-polystyrene` for MIE - a water suspension of polystyrene spheres.

    * `mie-fusedsilica` for MIE - a water suspension of fused silica spheres.

* `pf_param_1`: is the first parameter of the scattering phase function
  formatted with two decimal digits and using `_` as the decimal separator:

    * `g-<g>` for HG, e.g. `g-0_10` for *g*=0.1.

    * `g-<g>` for MHG, e.g. `g-0_50` for *g*=0.5.

    * `g-<g>` for GK, e.g. `g-0_90` for *g*=0.9.

    * `diameter-<d>um` MIE, e.g. `diameter-0_25` for *d*=0.25&nbsp;µm.

* `pf_param_2`: is the second parameter of the scattering phase function
  formatted with two decimal digits and using `_` as the decimal
  separator. An exception to this rule is the wavelength parameter of the
  MIE scattering phase function that is converted to nm and formatted as an
  integer. This placeholder is not used with the HG scattering phase function.

    * `b-<b>` for MHG, e.g. `b-0_60` for *β*=0.6.

    * `a-<a>` for GK, e.g. `a-0_50` for *α*=0.5.

    * `wavelength-<w>nm` for MIE, e.g. `wavelength-500nm` for *w*=500&nbsp;nm.

* `<mua>` is the absorption coefficient in units of cm<sup>-1</sup> with
  two decimal digits and `_` as the decimal separator, e.g `2_50`
  for *µ<sub>a</sub>*=2.5&nbsp;<sup>-1</sup>.

* `<musr>` is the reduced scattering coefficient in units of cm<sup>-1</sup>
  with two decimal digits and `_` as a decimal separator, e.g `20_00`
  for *µ<sub>a</sub>*=20.0&nbsp;cm<sup>-1</sup>.

## Voxelized media datasets (MCVOX)

These datasets include fluence/energy deposition data simulated with the
MC kernel for voxelized media. A two-layer skin model with an embedded blood
vessel is used. The depth/position of the blood vessel along the *z* axis is
varied from 0.2 to 0.8&nbsp;mm in steps of 0.025&nbsp;mm.
The refractive index of the surrounding medium is set to 1.337.
The simulations are run with 1000 million photon packets.

*A two-layer skin model with an embedded blood vessel*
|              | Parameter         | Value                                       | Description                        |
|--------------|-------------------|---------------------------------------------|------------------------------------|
| Line         |                   |                                             | normally incident                  |
| Material     | *µ<sub>a</sub>*   | 16.5724&nbsp;cm<sup>-1</sup>                | epidermis                          |
|              | *µ<sub>s</sub>*   | 375.9398&nbsp;cm<sup>-1</sup>               |                                    |
|              | *n*               | 1.337                                       |                                    |
|              | *pf*              | HG(0.9)                                     |                                    |
| Material     | *µ<sub>a</sub>*   | 45.85&nbsp;cm<sup>-1</sup>                  | dermis                             |
|              | *µ<sub>s</sub>*   | 356.5406&nbsp;cm<sup>-1</sup>               |                                    |
|              | *n*               | 1.337                                       |                                    |
|              | *pf*              | HG(0.9)                                     |                                    |
| Material     | *µ<sub>a</sub>*   | 230.5427&nbsp;cm<sup>-1</sup>               | blood vessel                       |
|              | *µ<sub>s</sub>*   | 93.985&nbsp;cm<sup>-1</sup>                 |                                    |
|              | *n*               | 1.337                                       |                                    |
|              | *pf*              | HG(0.9)                                     |                                    |
| Fluence      | *xaxis*           | Axis(start=-502.5e-6, stop=502.5e-6, n=201) | energy deposition detector         |
|              | *yaxis*           | Axis(start=-502.5e-6, stop=502.5e-6, n=201) |                                    |
|              | *zaxis*           | Axis(start=0.0, stop=0.001, n=200)          |                                    |
| Cartesian    | *xaxis*           | Axis(start=-502.5e-6, stop=502.5e-6, n=201) | reflectance/transmittance detector |
|              | *yaxis*           | Axis(start=-502.5e-6, stop=502.5e-6, n=201) |                                    |
| Blood vessel | *diameter*        | 200.0&nbsp;μm                               | in dermis                          |
|              | *position*        | (x, y) = (0, 0)                             |                                    |
|              | *z*               | from 0.2 to 0.8 in 0.025&nbsp;mm steps      |                                    | 
|              | *direction*       | (x, y, z) = (0, 1, 0)                       |                                    |
| Epidermis    | *thickness*       | 100&nbsp;μm                                 |                                    |

![Cartesian](/docs/source/static/mcdetector_cartesian.svg)<br>
*Illustration of the Cartesian detector that is used to collect the reflectance and transmittance.*

![Deposition](/docs/source/static/deposition-projection.gif)<br>
*Deposition mean along the y axis for the above configuration of a 2-layer skin with an embedded blood vessel*

### Dataset files

Datasets are available as compressed numpy data files that are organized as
follows:

```
data/mcvox/fluence/2-layer-skin-<diameter>um-vessel-<depth>um-depth-deposition.npz
```

The values of placeholders <> are as follows:

* `<diameter>` is the diameter of the blood vessel in units of μm, formatted
  as an integer value, e.g `200` for a 200&nbsp;μm blood vessel.

* `<depth>` is the *z* coordinate (depth) of the blood vessel in units
  of μm, formatted as an integer value, e.g `500` for *z*=500&nbsp;μm.

## Sampling volume datasets (SV)

The sampling volume dataset is computed for a semi-infinite homogeneous medium
for an optical fiber probe with two optical fibers placed at a distance of
0.5&nbsp;mm. The refractive index of the surrounding medium is set to 1.0.
Simulations are run in batches until 1,000,000 photon packet traces that reach
the detector fiber are collected and converted to sampling volume information.
The trace capacity is limited to 1000 events. The simulation termination
radius is set to 25&nbsp;mm. 

*Sampling volume for a probe with two optical fibers*
|                | Parameter        | Value                                       | Description                  |
|----------------|------------------|---------------------------------------------|------------------------------|
| LinearArray    | *dcore*          | 200&nbsp;μm                                 | linear layout of 2 fibers    |
|                | *dcladding*      | 220&nbsp;μm                                 |                              |
|                | *ncore*          | 1.462                                       |                              |
|                | *NA*             | 0.22                                        |                              |
|                | *spacing*        | 500&nbsp;μm                                 |                              |
|                | *diameter*       | 6.0&mm                                      |                              |
|                | *reflectivity*   | 0.6                                         |                              |
| Trace          | *maxlen*         | 1000                                        | packet trace configuration   |
| SamplingVolume | *xaxis*          | Axis(start=-0.00075, stop=0.00075, n=300)   |                              |
|                | *yaxis*          | Axis(start=-0.00075, stop=0.00075, n=300)   |                              |
|                | *zaxis*          | Axis(start=0.0, stop=0.001, n=200)          | sampling volume voxelization |
| Layer          | *μ<sub>a</sub>*  | 2.0&nbsp;cm<sup>-1</sup>                    | sample layer                 |
|                | *μ<sub>s</sub>*  | 500&nbsp;cm<sup>-1</sup>                    |                              |
|                | *n*              | 1.337                                       |                              |
|                | *pf*             | HG(0.95)                                    |                              |

![Sampling volume](/docs/source/static/sv-200um-0_22na-500um-sds-projection.svg)<br>
*Sampling volume mean along the y axis for the above configuration of optical probe with two optical fibers*


### Dataset files

Datasets are available as compressed numpy data files that are organized as
follows:

```
data/mcml/1-layer-semiinfinite/sv/reflectance/fiber-200um-0_22na/sds-<sds>um/hg/g-<g>/um-mua-<mua>-musr-<musr>-invcm.npz
```

The values of placeholders <> are as follows:

* `<sds>` is the distance between the centers of the source and detector
  fibers in units of μm, formatted as an integer value, e.g `500` for a 
  500&nbsp;μm distance.

* `<g>` is the anisotropy formatted with two decimal digits and `_`
  as the decimal separator, e.g `0_15` for *g*=0.15.

* `<mua>` is the absorption coefficient in units of cm<sup>-1</sup> with
  two decimal digits and `_` as the decimal separator, e.g `2_50`
  for *μ<sub>a</sub>*=2.5&nbsp;cm<sup>-1</sup>.

* `<musr>` is the reduced scattering coefficient in units of cm<sup>-1</sup>
  with two decimal digits and `_` as a decimal separator, e.g `20_00`
  for *μ<sub>s</sub>'*=20.0&nbsp;cm<sup>-1</sup>.

## Citing MC Dataset

We, the authors of MC Dataset, expect that the package is used in accordance
with the [GPL3+](https://www.gnu.org/licenses/gpl-3.0-standalone.html)
license and that any work using the MC Dataset package also cites the project
and at least one of the following references:

 - P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Limitations of the commonly*
   *used simplified laterally uniform optical fiber probe-tissue interface in*
   *Monte Carlo simulations of diffuse reflectance*, Biomed. Opt. Expres,
   **6** (10), 3973-3988 (2015), https://doi.org/10.1364/BOE.6.003973.

 - P. Naglič, F. Pernuš, B. Likar, and M. Bürmen, *Lookup table-based sampling*
   *of the phase function for Monte Carlo simulations of light propagation in*
   *turbid media*, Biomed. Opt. Expres, **8** (3), 1895-1910 (2017),
   https://doi.org/10.1364/BOE.8.001895.

For alternative licensing options of MC Dataset please contact us at
info@xopto.eu.
