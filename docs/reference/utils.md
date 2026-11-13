Utility functions and classes.

## Parallelism

The [Queue][dexter.Queue] type can utilize multiple threads when running particle routines. To set the number of threads use [`set_num_threads`][dexter.utils.set_num_threads]:

::: dexter.utils.set_num_threads
    options:
      heading_level: 3

To get the number of threads available on your system, use [`get_max_threads`][dexter.utils.get_max_threads]:

::: dexter.utils.get_max_threads
    options:
      heading_level: 3

## Colors

::: dexter.colors.particle_color
    options:
      heading_level: 3
      show_source: true

::: dexter.colors.orbit_color
    options:
      heading_level: 3
      show_source: true

::: dexter.colors.critical_point_color
    options:
      heading_level: 3
      show_source: true

::: dexter.colors.boundary_extremum_color
    options:
      heading_level: 3
      show_source: true
