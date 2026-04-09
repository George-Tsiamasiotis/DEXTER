# Parallelism

Using Rust's [rayon] crate, particle routines can run in parallel. The number of threads to be used can be adjusted with the [`set_num_threads`][dexter.set_num_threads] function.

By default, Rust will use all the available threads. If this results in the device being unresponsive, consider lowering the max thread number. The device's number of available threads can be found with [`get_max_threads`][dexter.get_max_threads].

::: dexter.set_num_threads

::: dexter.get_max_threads

[rayon]: https://docs.rs/rayon/latest/rayon/
