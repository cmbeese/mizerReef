![mizerReef logo](reference/figures/logo.png)

## Overview

MizerReef is an R package that extends the
[mizer](https://github.com/sizespectrum/mizer) framework to support
modeling of structurally complex marine ecosystems such as coral reefs,
kelp forests, and rocky reefs. Predator-prey interactions are mediated
by predation refuge to explicitly capture benthic habitat structure. It
also adds features for benthic coupling, introducing additional
resources such as detritus and algae that are not size-structured and
can be consumed by all relevant species regardless of size. MizerReef is
suitable for researchers studying any marine habitat where structural
complexity influences ecosystem dynamics.

> **Compatibility note:** MizerReef is not yet verified to work
> correctly when loaded together with other mizer extension packages
> (e.g. mizerMR). A mechanism meant to keep bundled example models
> correctly classed in that situation is not currently active; this is
> being worked on. See
> [`vignette("extension_mechanism")`](https://cmbeese.github.io/mizerReef/articles/extension_mechanism.md)
> for details. MizerReef on its own is unaffected.

MizerReef was originally developed to support the creation and
exploration of a generic model for tropical coral reefs as part of a PhD
thesis. Example parameters for this system are included with the
package. For details, see: [Modelling Coral Reef Futures: Exploring the
role of structural complexity in sustaining ecosystem
services](https://doi.org/10.26686/wgtn.26421523).

## Installation

MizerReef is under active development. Feedback is welcome! You can
install the development version of MizerReef from GitHub:

``` r

# If you don't have the remotes package, install it first:
install.packages("remotes")

# Then install MizerReef:
remotes::install_github("cmbeese/mizerReef")
```

Make sure you have the latest versions of
[mizer](https://github.com/sizespectrum/mizer) and
[mizerExperimental](https://github.com/sizespectrum/mizerExperimental),
as MizerReef is designed to work with their most recent releases.

## Vignettes and Examples

See the [mizerReef website](https://cmbeese.github.io/mizerReef/) for
tutorials and example workflows.

## Getting Help

For questions, bug reports, or support, please open an issue on
[GitHub](https://github.com/cmbeese/mizerReef/issues).

## Citation

If you use MizerReef in your research, please cite:

Beese C. M., Delius G., Mumby P. J., Rogers A. (2026). *MizerReef: an R
package to create multi-species size spectrum models for structurally
complex habitats*. R package version 2.0.0,
<https://github.com/cmbeese/mizerReef>.

> A peer-reviewed publication and DOI for MizerReef will be available
> soon. Please check back for the updated citation.

## License

MizerReef is released under the GPL-3 license.
