//! Methods related to binary search through a slice of floats.

/// Returns the index i of the array `xarray` such that `xarray[i] <= x <= xarray[i+1]`.
///
/// The index is searched for in the range `[ilo, ihi]`, which can be used to speed up the search.
///
/// Code was taken from [`rsl_interpolation's`] [`Accelerator's`][`rsl_interpolation::Accelerator`]
/// `bsearch` method.
///
/// [`rsl_interpolation's`]: https://github.com/George-Tsiamasiotis/rsl-interpolation
#[must_use]
#[expect(clippy::shadow_same, reason = "clarity, gsl")]
pub fn bsearch(xarray: &[f64], x: f64, ilo: usize, ihi: usize) -> usize {
    let mut ilo = ilo;
    let mut ihi = ihi;
    while ihi > ilo + 1 {
        let mid = ihi.midpoint(ilo);
        if xarray[mid] > x {
            ihi = mid;
        } else {
            ilo = mid;
        }
    }
    ilo
}

/// Performs a lookup action on a data array. Returns an index `i` such that
/// `xarray[i] <= x < xarray[i+1]`.
///
/// The `cache` parameter is used as an initial guess to speed up the [`binary search`][bsearch].
///
/// # Panics
///
/// This method panics if `cache` is larger that `xa`'s length.
#[must_use]
pub fn find(xarray: &[f64], x: f64, cache: usize) -> usize {
    if x < xarray[cache] {
        bsearch(xarray, x, 0, cache)
    } else if x >= xarray[cache + 1] {
        bsearch(xarray, x, cache, xarray.len() - 1)
    } else {
        cache
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn setup_xarray() -> [f64; 5] {
        [0.0, 1.0, 2.0, 3.0, 4.0]
    }

    #[test]
    fn bsearch_interior_point() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, 1.5, 0, 4);
        assert_eq!(res, 1);
    }

    #[test]
    fn search_last_value() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, 4.0, 0, 4);
        assert_eq!(res, 3);
    }

    #[test]
    fn search_first_value() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, 0.0, 0, 4);
        assert_eq!(res, 0);
    }

    #[test]
    fn bsearch_boundary() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, 2.0, 0, 4);
        assert_eq!(res, 2);
    }

    #[test]
    fn bsearch_above_bounds() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, 10.0, 0, 4);
        assert_eq!(res, 3);
    }

    #[test]
    fn bsearch_below_bounds() {
        let xarray = setup_xarray();

        let res = bsearch(&xarray, -10.0, 0, 4);
        assert_eq!(res, 0);
    }

    #[test]
    fn bsearch_find() {
        let xarray = setup_xarray();
        let mut k1 = 0;
        let mut k2 = 0;
        let mut t = false;
        let mut cache = 0;
        let r = [
            -0.2, 0.0, 0.1, 0.7, 1.0, 1.3, 1.9, 2.0, 2.2, 2.7, 3.0, 3.1, 3.6, 4.0, 4.1, 4.9,
        ];

        // run through all the pairs of points.
        while (k1 < 16) & (k2 < 16) {
            let x = if t { r[k1] } else { r[k2] };
            t = !t;

            if !t {
                k1 = (k1 + 1) % 16;
                if k1 == 0 {
                    k2 += 1;
                }
            }

            cache = find(&xarray, x, cache);
            let j = bsearch(&xarray, x, 0, 4);
            assert_eq!(cache, j);
        }
    }
}
