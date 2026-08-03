# Get energy rate available for growth through time

Calculates the energy rate \\g_i(w)\\ (grams/year) available by species
and size for growth after metabolism, movement and reproduction have
been accounted for.

## Usage

``` r
getEGrowthTime(object, n, n_pp, n_other, time_range, drop = FALSE, ...)
```

## Arguments

- object:

  A `MizerParams` object or a `MizerSim` object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- time_range:

  A numeric or character vector of times. Only the range of values
  matters, so all saved times between `min(time_range)` and
  `max(time_range)` are selected.

- drop:

  If `drop = TRUE` then the dimension of length 1 will be removed from
  the returned array.

- ...:

  Unused

## Value

If a `MizerParams` object is passed in, the function returns a two
dimensional array (predator species x predator size) based on the
abundances also passed in. If a `MizerSim` object is passed in, the
function returns a three dimensional array (time step x predator species
x predator size) with the energy for growth calculated at every time
step in the simulation. If `drop = TRUE` then the dimension of length 1
will be removed from the returned array.

## See also

[`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md)
