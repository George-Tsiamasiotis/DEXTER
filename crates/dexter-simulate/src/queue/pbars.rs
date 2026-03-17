//! Progress bar customization for Queue's routines.

#![expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering::SeqCst};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

use crate::Queue;
use crate::{IntegrationStatus, IntersectParams};

/// The [`Queue::integrate`] progress bar style.
const INTEGRATE_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.cyan/blue}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The [`Queue::integrate`] progress bar chars (filled, current, to do).
const INTEGRATE_PROGRESS_CHARS: &str = "#>-";

/// The [`Queue::intersect`] progress bar style.
const INTERSECT_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.cyan/blue}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The [`Queue::intersect`] progress bar chars (filled, current, to do).
const INTERSECT_PROGRESS_CHARS: &str = "#>-";

// ===============================================================================================
// ===============================================================================================

/// [`Queue::integrate`] progress bar helper struct.
pub(crate) struct IntegratePbar {
    pbar: ProgressBar,
    length: usize,
    // Live statistics
    out_of_bounds: Arc<AtomicUsize>,
    integrated: Arc<AtomicUsize>,
    escaped: Arc<AtomicUsize>,
    timed_out: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
}

impl IntegratePbar {
    /// Initialize the pre-configured progress bar.
    pub(crate) fn new(queue: &Queue) -> Self {
        let style = ProgressStyle::with_template(INTEGRATE_PBAR_STYLE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars(INTEGRATE_PROGRESS_CHARS);
        let pbar = ProgressBar::new(queue.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        Self {
            pbar,
            length: queue.particles.len(),
            out_of_bounds: Arc::default(),
            integrated: Arc::default(),
            escaped: Arc::default(),
            timed_out: Arc::default(),
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
        self.pbar.println("🗿 Integrating");
    }

    /// Increases the wrapped pbar, as well as the live statistics.
    pub(crate) fn inc(&self, status: &IntegrationStatus) {
        self.pbar.inc(1);
        let _: usize = match *status {
            IntegrationStatus::OutOfBoundsInitialization => self.out_of_bounds.fetch_add(1, SeqCst),
            IntegrationStatus::Integrated => self.integrated.fetch_add(1, SeqCst),
            IntegrationStatus::Escaped => self.escaped.fetch_add(1, SeqCst),
            IntegrationStatus::TimedOut(..) => self.timed_out.fetch_add(1, SeqCst),
            IntegrationStatus::Failed(..) => self.failed.fetch_add(1, SeqCst),
            _ => unreachable!("no other possible status"),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "===== 📊 Stats =====\n",
                "ℹ️ OutOfBounds = {}\n",
                "✅ Integrated  = {}\n",
                "🏃 Escaped     = {}\n",
                "⌛ Timed-out   = {}\n",
                "🥀 Failed      = {}",
            ),
            self.out_of_bounds.load(SeqCst),
            self.integrated.load(SeqCst),
            self.escaped.load(SeqCst),
            self.timed_out.load(SeqCst),
            self.failed.load(SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("👌 Integration Done");
        self.pbar.force_draw();
        self.pbar.finish();
    }
}

// ===============================================================================================

/// [`Queue::intersect`] progress bar helper struct.
pub(crate) struct IntersectPbar {
    pbar: ProgressBar,
    length: usize,
    intersect_params: IntersectParams,
    // Live statistics
    out_of_bounds: Arc<AtomicUsize>,
    intersected: Arc<AtomicUsize>,
    intersected_timed_out: Arc<AtomicUsize>,
    escaped: Arc<AtomicUsize>,
    timed_out: Arc<AtomicUsize>,
    invalid_intersections: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
}

impl IntersectPbar {
    /// Initialize the pre-configured progress bar.
    pub(crate) fn new(queue: &Queue, intersect_params: &IntersectParams) -> Self {
        let style = ProgressStyle::with_template(INTERSECT_PBAR_STYLE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars(INTERSECT_PROGRESS_CHARS);
        let pbar = ProgressBar::new(queue.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        Self {
            pbar,
            length: queue.particles.len(),
            intersect_params: intersect_params.clone(),
            out_of_bounds: Arc::default(),
            intersected: Arc::default(),
            intersected_timed_out: Arc::default(),
            escaped: Arc::default(),
            timed_out: Arc::default(),
            invalid_intersections: Arc::default(),
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
            "🗿 Integrating with {:?}={:.4} for {} turns",
            self.intersect_params.intersection,
            self.intersect_params.angle,
            self.intersect_params.turns,
        ));
    }

    /// Increases the wrapped pbar, as well as the live statistics.
    pub(crate) fn inc(&self, status: &IntegrationStatus) {
        self.pbar.inc(1);
        #[rustfmt::skip]
        let _: usize = match *status {
            IntegrationStatus::OutOfBoundsInitialization => self.out_of_bounds.fetch_add(1, SeqCst),
            IntegrationStatus::Intersected => self.intersected.fetch_add(1, SeqCst),
            IntegrationStatus::IntersectedTimedOut => self.intersected_timed_out.fetch_add(1, SeqCst),
            IntegrationStatus::Escaped => self.escaped.fetch_add(1, SeqCst),
            IntegrationStatus::TimedOut(..) => self.timed_out.fetch_add(1, SeqCst),
            IntegrationStatus::InvalidIntersections => self.invalid_intersections.fetch_add(1, SeqCst),
            IntegrationStatus::Failed(..) => self.failed.fetch_add(1, SeqCst),
            _ => unreachable!("no other possible status"),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "========= 📊 Stats =========\n",
                "ℹ️ OutOfBounds         = {}\n",
                "✅ Intersected         = {}\n",
                "⌛ IntersectedTimedOut = {}\n",
                "🏃 Escaped             = {}\n",
                "⌛ Timed-out           = {}\n",
                "⛔ Invalid             = {}\n",
                "🥀 Failed              = {}",
            ),
            self.out_of_bounds.load(SeqCst),
            self.intersected.load(SeqCst),
            self.intersected_timed_out.load(SeqCst),
            self.escaped.load(SeqCst),
            self.timed_out.load(SeqCst),
            self.invalid_intersections.load(SeqCst),
            self.failed.load(SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("👌 Intersection Done");
        self.pbar.force_draw();
        self.pbar.finish();
    }
}
