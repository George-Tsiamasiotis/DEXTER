//! Progress bar styles and methods for Heap calculations.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};
use particle::{IntegrationStatus, IntegrationStatus::*, MappingParameters};

use crate::Heap;

/// The Poincare map calculation progress bar style.
pub const POINCARE_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.cyan/blue}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The Poincare map progress bar chars (filled, current, to do).
pub const POINCARE_PROGRESS_CHARS: &str = "#>-";

/// The frequency calculation progress bar style.
pub const FREQUENCY_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.orange/red}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The Poincare map progress bar chars (filled, current, to do).
pub const FREQUENCY_PROGRESS_CHARS: &str = "#>-";

// ===============================================================================================

pub(crate) struct PoincarePbar {
    pbar: ProgressBar,
    params: MappingParameters,
    length: usize,
    // Live statistics
    mapped: Arc<AtomicUsize>,
    escaped: Arc<AtomicUsize>,
    timedout: Arc<AtomicUsize>,
    invalid: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
}

impl PoincarePbar {
    /// Initializes the progress bar.
    pub(crate) fn new(heap: &Heap, params: &MappingParameters) -> Self {
        // `.progress_with()` seems to update the pbar **before** the map() method is called, so
        // we must create it and update it manually.
        let style = ProgressStyle::with_template(POINCARE_PBAR_STYLE)
            .unwrap_or(ProgressStyle::default_bar())
            .progress_chars(POINCARE_PROGRESS_CHARS);
        let pbar = ProgressBar::new(heap.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        Self {
            pbar,
            params: *params,
            length: heap.particles.len(),
            mapped: Arc::default(),
            escaped: Arc::default(),
            timedout: Arc::default(),
            invalid: Arc::default(),
            failed: Arc::default(),
        }
    }

    /// Prints an informative message before the ticking starts.
    pub(crate) fn print_prelude(&self) {
        self.pbar.println(format!(
            "🚀 Using {} threads for {} particles",
            rayon::current_num_threads(),
            self.length
        ));
        self.pbar.println(format!(
            "🗿 Integrating with {:?}={:.4} for {} intersections",
            self.params.section, self.params.alpha, self.params.intersections,
        ));
    }

    /// Increases the wrapped pbar, as well as the live statistics
    pub(crate) fn inc(&self, status: &IntegrationStatus) {
        self.pbar.inc(1);
        match status {
            Mapped => self.mapped.fetch_add(1, Ordering::SeqCst),
            Escaped => self.escaped.fetch_add(1, Ordering::SeqCst),
            TimedOut(..) => self.timedout.fetch_add(1, Ordering::SeqCst),
            InvalidIntersections => self.invalid.fetch_add(1, Ordering::SeqCst),
            Failed(..) => self.failed.fetch_add(1, Ordering::SeqCst),
            _ => 0, // ignored
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                // "📊 Stats:\n",
                "📍 Mapped    = {}\n",
                "🏃 Escaped   = {}\n",
                "⌛ Timed-out = {}\n",
                "⛔ Invalid   = {}\n",
                "🥀 Failed    = {}",
            ),
            self.mapped.load(Ordering::SeqCst),
            self.escaped.load(Ordering::SeqCst),
            self.timedout.load(Ordering::SeqCst),
            self.invalid.load(Ordering::SeqCst),
            self.failed.load(Ordering::SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("✅️ Integration Done");
        self.pbar.finish();
    }
}

// ===============================================================================================

pub(crate) struct FrequenciesPbar {
    pbar: ProgressBar,
    length: usize,
    // Live statistics
    integrated: Arc<AtomicUsize>,
    escaped: Arc<AtomicUsize>,
    timedout: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
}

impl FrequenciesPbar {
    /// Initializes the progress bar.
    pub(crate) fn new(heap: &Heap) -> Self {
        // `.progress_with()` seems to update the pbar **before** the map() method is called, so
        // we must create it and update it manually.
        let style = ProgressStyle::with_template(FREQUENCY_PBAR_STYLE)
            .unwrap_or(ProgressStyle::default_bar())
            .progress_chars(FREQUENCY_PROGRESS_CHARS);
        let pbar = ProgressBar::new(heap.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        Self {
            pbar,
            length: heap.particles.len(),
            integrated: Arc::default(),
            escaped: Arc::default(),
            timedout: Arc::default(),
            failed: Arc::default(),
        }
    }

    /// Prints an informative message before the ticking starts.
    pub(crate) fn print_prelude(&self) {
        self.pbar.println(format!(
            "🚀 Using {} threads for {} particles",
            rayon::current_num_threads(),
            self.length
        ));
    }

    /// Increases the wrapped pbar, as well as the live statistics
    pub(crate) fn inc(&self, status: &IntegrationStatus) {
        self.pbar.inc(1);
        match status {
            SinglePeriodIntegrated => self.integrated.fetch_add(1, Ordering::SeqCst),
            Escaped => self.escaped.fetch_add(1, Ordering::SeqCst),
            TimedOut(..) => self.timedout.fetch_add(1, Ordering::SeqCst),
            Failed(..) => self.failed.fetch_add(1, Ordering::SeqCst),
            _ => 0, // ignored
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                // "📊 Stats:\n",
                "👌 Closed    = {}\n",
                "🏃 Escaped   = {}\n",
                "⌛ Timed-out = {}\n",
                "🥀 Failed    = {}",
            ),
            self.integrated.load(Ordering::SeqCst),
            self.escaped.load(Ordering::SeqCst),
            self.timedout.load(Ordering::SeqCst),
            self.failed.load(Ordering::SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("✅️ Integration Done");
        self.pbar.finish();
    }
}
