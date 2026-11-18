# Changelog

## [0.2.8] - DATE

### Added

- Add electric potential evaluation methods for `Mode` trait and `Perturbation` type.
- Add `dΦ/dθ` and `dΦ/dζ` terms in the GC equations of motion.

### Changed

- Refactor `Mode` and `Perturbation` evaluation methods' names to include `m` and `phi`.
- Refactor perturbation-related terms on a `PerturbationTerms` struct in `GCState`.
- Split the electric potential on the integrable and perturbative components in `GCState`.

[0.2.8]: https://github.com/George-Tsiamasiotis/DEXTER/releases/tag/v0.2.8
