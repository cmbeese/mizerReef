# Contribution of unstructured components to the encounter rate

The encounter rate \\E_i(w)\\ for an unstructured resource like algae or
detritus is proportional to the total biomass \\B\\ with a coefficient
\\\rho_i(w)\\ that depends on the species \\i\\ and the size of the
consumer:

## Usage

``` r
encounter_contribution(params, n_other, component, ...)
```

## Arguments

- params:

  MizerParams

- n_other:

  Biomasses of unstructured components

- component:

  Name of component whose contribution is requested

- ...:

  Unused

## Value

Array (species x size) with the encounter rate in g/year.

## Details

\$\$E_i(w) = \rho_i(w) B.\$\$

The coefficient \\\rho_i(w)\\ is stored as a matrix (species x size) in
the `rho` parameter of the component. It has units 1/year.
