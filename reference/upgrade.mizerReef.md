# Upgrade a mizerReef params object to the current layout

Called automatically by
[`mizer::validParams()`](https://sizespectrum.org/mizer/reference/validParams.html)
when an object created with an older version of mizerReef is loaded.
Migrates data stored in the old custom S4 slots (`refuge_params`,
`algae_params`, `detritus_params`) to the current `other_params`
sub-lists (`other_params$refuge_params`, `other_params$algae`,
`other_params$detritus`).

## Usage

``` r
# S3 method for class 'mizerReef'
upgrade(object, ...)
```

## Arguments

- object:

  A `mizerReef` params object (possibly from an older version).

- ...:

  Unused.

## Value

An upgraded `mizerReef` object.
