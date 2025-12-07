//! Representation of a total perturbation, a sum of multiple harmonics.

use rsl_interpolation::Accelerator;

use crate::cache::HarmonicCache;
use crate::harmonics::NcHarmonic;
use crate::{Harmonic, Perturbation};

/// A sum of different perturbation harmonics.
pub struct NcPerturbation {
    harmonics: Vec<NcHarmonic>,
}

// Creation and data extraction
impl NcPerturbation {
    pub fn from_harmonics(harmonics: &[NcHarmonic]) -> Self {
        Self {
            harmonics: harmonics.into(),
        }
    }

    pub fn get_harmonics(&self) -> Vec<NcHarmonic> {
        self.harmonics.clone()
    }
}

impl Perturbation for NcPerturbation {
    fn p(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .h(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dpsip(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dtheta(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dzeta(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dt(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }
}

/// Getters
impl NcPerturbation {
    /// Returns the number of harmonics.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.harmonics.len()
    }
}

impl std::fmt::Debug for NcPerturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for harmonic in self.get_harmonics() {
            let _ = harmonic.fmt(f);
        }
        Ok(())
    }
}

// #[cfg(test)]
// mod test {
//     use crate::extract::STUB_TEST_NETCDF_PATH;
//     use crate::harmonics::NcHarmonicBuilder;
//     use crate::*;
//     use rsl_interpolation::Accelerator;
//
//     use std::path::PathBuf;
//
//     #[test]
//     fn test_summation() {
//         let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
//         let per1 = NcPerturbation::from_harmonics(&vec![
//             NcHarmonicBuilder::new(&path, "akima", 2, 1)
//                 .build()
//                 .unwrap(),
//         ]);
//         let per2 = NcPerturbation::from_harmonics(&vec![
//             NcHarmonicBuilder::new(&path, "akima", 3, 1)
//                 .build()
//                 .unwrap(),
//             NcHarmonicBuilder::new(&path, "akima", 3, 2)
//                 .build()
//                 .unwrap(),
//             NcHarmonicBuilder::new(&path, "akima", 2, 1)
//                 .build()
//                 .unwrap(),
//         ]);
//
//         let mut acc = Accelerator::new();
//         // Normally, this would happen inside State.
//         let mut hcache1: Vec<Box<dyn HarmonicCache>> = vec![Box::new(HarmonicCache::default())];
//         let mut hcache2: Vec<Box<dyn HarmonicCache>> = vec![Box::new(HarmonicCache::default())];
//         let psip = 0.000001;
//         let theta = 1.0;
//         let zeta = 1.0;
//
//         assert_eq!(
//             3.0 * per1.p(psip, theta, zeta, &mut hcache1, &mut acc).unwrap(),
//             per2.p(psip, theta, zeta, &mut hcache2, &mut acc).unwrap(),
//         );
//         assert_eq!(
//             3.0 * per1
//                 .dp_dpsip(psip, theta, zeta, &mut hcache1, &mut acc)
//                 .unwrap(),
//             per2.dp_dpsip(psip, theta, zeta, &mut hcache2, &mut acc)
//                 .unwrap(),
//         );
//         assert_eq!(
//             3.0 * per1
//                 .dp_dtheta(psip, theta, zeta, &mut hcache1, &mut acc)
//                 .unwrap(),
//             per2.dp_dtheta(psip, theta, zeta, &mut hcache2, &mut acc)
//                 .unwrap(),
//         );
//         assert_eq!(
//             3.0 * per1
//                 .dp_dzeta(psip, theta, zeta, &mut hcache1, &mut acc)
//                 .unwrap(),
//             per2.dp_dzeta(psip, theta, zeta, &mut hcache2, &mut acc)
//                 .unwrap(),
//         );
//         assert_eq!(
//             3.0 * per1
//                 .dp_dt(psip, theta, zeta, &mut hcache1, &mut acc)
//                 .unwrap(),
//             per2.dp_dt(psip, theta, zeta, &mut hcache2, &mut acc)
//                 .unwrap(),
//         );
//     }
//
//     #[test]
//     fn test_perturbation_misc() {
//         let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
//         let harmonics = vec![
//             NcHarmonicBuilder::new(&path, "akima", 2, 1)
//                 .build()
//                 .unwrap(),
//         ];
//         let per = NcPerturbation::from_harmonics(&harmonics);
//         assert_eq!(per.len(), 1);
//
//         let _ = per.get_harmonics();
//         let _ = format!("{:?}", per);
//     }
// }
