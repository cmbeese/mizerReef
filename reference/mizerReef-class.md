# mizerReef extension classes

S3 extension classes for
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
and [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
that enable S3 dispatch for extension-specific methods.

## Details

The class names are ordinary entries in the object's S3 class vector.
All reef-specific data lives in `other_params(params)` or in component
parameters (see
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)).

Objects of class `mizerReef` are created by
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md).
Objects of class `mizerReefSim` are returned automatically by
[`project()`](https://sizespectrum.org/mizer/reference/project.html)
when called on a `mizerReef` params object.

No class declaration is needed.
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
records the extension on the object with
[`mizer::recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.html)
and then calls
[`mizer::coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.html).

## See also

[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html),
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
