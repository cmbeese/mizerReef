# S4 marker class for mizerReef extension

A formal S4 class that extends
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
but adds no new slots. All reef-specific state (refuge, algae, detritus
parameters) lives in `other_params(params)`. The class label enables S3
dispatch of mizer generics to mizerReef-specific methods.

## Details

The class is **not** defined statically. Instead mizer creates it when
the package is loaded: `.onLoad()` calls
[`mizer::registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.html),
which recognises mizerReef as a dispatching extension from the S3
methods it registers for its marker class and inserts `mizerReef` at the
correct place in the S4 hierarchy relative to any other extension
packages loaded in the same session. This lets mizerReef be chained with
other extensions (for example mizerMR, to add multiple resources) in
either load order. A static `contains = "MizerParams"` definition would
fix mizerReef as a direct sibling of every other extension and prevent
such chaining, because a sealed class cannot be re-parented.

## See also

[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html),
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
