//! Progress bar customization for Queue's routines.

#![expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]

use colored::Colorize;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering::SeqCst};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

use dexter_comspace::OrbitType;

use crate::{IntegrationStatus, IntersectParams, Queue};

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

/// The [`Queue::close`] progress bar style.
const CLOSE_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.cyan/blue}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The [`Queue::close`] progress bar chars (filled, current, to do).
const CLOSE_PROGRESS_CHARS: &str = "#>-";

/// The [`Queue::classify`] progress bar style.
const CLASSIFY_PBAR_STYLE: &str = concat!(
    "{msg}\n", // for Stats
    "🕜 {elapsed_precise} ",
    "{prefix} ",
    "[{wide_bar:.cyan/blue}] ",
    "{spinner:.bold} ",
    "{pos:>2}/{len:2} ",
    "({eta}) ",
);

/// The [`Queue::close`] progress bar chars (filled, current, to do).
const CLASSIFY_PROGRESS_CHARS: &str = "#>-";

mod orbit_type_colors {
    pub(super) const UNDEFINED: &str = "#fc5a50";
    pub(super) const TRAPPED_CONFINED: &str = "#ff000d";
    pub(super) const TRAPPED_LOST: &str = "#9a0200";
    pub(super) const COPASSING_CONFINED: &str = "#0165fc";
    pub(super) const COPASSING_LOST: &str = "#040273";
    pub(super) const CUPASSING_CONFINED: &str = "#01ff07";
    pub(super) const CUPASSING_LOST: &str = "#02590f";
    pub(super) const POTATO: &str = "#d1b26f";
    pub(super) const STAGNATED: &str = "#82cafc";
    pub(super) const UNCLASSIFIED: &str = "#be03fd";
}

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
    green: String,
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
            green: "Integrated".green().bold().to_string(),
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
            _ => self.failed.fetch_add(1, SeqCst),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "===== 📊 Stats =====\n",
                "✅ {}  = {}\n",
                "🧱 OutOfBounds = {}\n",
                "🏃 Escaped     = {}\n",
                "⌛ Timed-out   = {}\n",
                "🥀 Failed      = {}",
            ),
            self.green,
            self.integrated.load(SeqCst),
            self.out_of_bounds.load(SeqCst),
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
    green: String,
}

impl IntersectPbar {
    /// Initialize the pre-configured progress bar.
    pub(crate) fn new(queue: &Queue, intersect_params: &IntersectParams) -> Self {
        let style = ProgressStyle::with_template(INTERSECT_PBAR_STYLE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars(INTERSECT_PROGRESS_CHARS);
        let pbar = ProgressBar::new(queue.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(200));
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
            green: "Intersected".green().bold().to_string(),
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
            _ => self.failed.fetch_add(1, SeqCst),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "========= 📊 Stats =========\n",
                "✅ {}         = {}\n",
                "🧱 OutOfBounds         = {}\n",
                "⌛ IntersectedTimedOut = {}\n",
                "🏃 Escaped             = {}\n",
                "⌛ Timed-out           = {}\n",
                "⛔ Invalid             = {}\n",
                "🥀 Failed              = {}",
            ),
            self.green,
            self.intersected.load(SeqCst),
            self.out_of_bounds.load(SeqCst),
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

// ===============================================================================================
// ===============================================================================================

/// [`Queue::close`] progress bar helper struct.
pub(crate) struct ClosePbar {
    pbar: ProgressBar,
    length: usize,
    // Live statistics
    out_of_bounds: Arc<AtomicUsize>,
    closed_periods: Arc<AtomicUsize>,
    escaped: Arc<AtomicUsize>,
    timed_out: Arc<AtomicUsize>,
    failed: Arc<AtomicUsize>,
    green: String,
}

impl ClosePbar {
    /// Initialize the pre-configured progress bar.
    pub(crate) fn new(queue: &Queue) -> Self {
        let style = ProgressStyle::with_template(CLOSE_PBAR_STYLE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars(CLOSE_PROGRESS_CHARS);
        let pbar = ProgressBar::new(queue.particles.len() as u64).with_style(style);
        Self {
            pbar,
            length: queue.particles.len(),
            out_of_bounds: Arc::default(),
            closed_periods: Arc::default(),
            escaped: Arc::default(),
            timed_out: Arc::default(),
            failed: Arc::default(),
            green: "Closed".green().bold().to_string(),
        }
    }

    /// Prints an informative message before the ticking starts.
    pub(crate) fn print_prelude(&self) {
        self.pbar.println(format!(
            "🚀 Using {} threads for {} particles",
            rayon::current_num_threads(),
            self.length
        ));
        self.pbar.println("🗿 Closing orbits");
    }

    /// Increases the wrapped pbar, as well as the live statistics.
    pub(crate) fn inc(&self, status: &IntegrationStatus) {
        self.pbar.inc(1);
        let _: usize = match *status {
            IntegrationStatus::OutOfBoundsInitialization => self.out_of_bounds.fetch_add(1, SeqCst),
            IntegrationStatus::ClosedPeriods(..) => self.closed_periods.fetch_add(1, SeqCst),
            IntegrationStatus::Escaped => self.escaped.fetch_add(1, SeqCst),
            IntegrationStatus::TimedOut(..) => self.timed_out.fetch_add(1, SeqCst),
            _ => self.failed.fetch_add(1, SeqCst),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "====== 📊 Stats =====\n",
                "✅ {}      = {}\n",
                "🧱 OutOfBounds = {}\n",
                "🏃 Escaped     = {}\n",
                "⌛ Timed-out   = {}\n",
                "🥀 Failed      = {}",
            ),
            self.green,
            self.closed_periods.load(SeqCst),
            self.out_of_bounds.load(SeqCst),
            self.escaped.load(SeqCst),
            self.timed_out.load(SeqCst),
            self.failed.load(SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("👌 Period closing done");
        self.pbar.force_draw();
        self.pbar.finish();
    }
}

// ===============================================================================================
// ===============================================================================================

/// [`Queue::classify`] progress bar helper struct.
pub(crate) struct ClassifyPbar {
    pbar: ProgressBar,
    length: usize,
    // Live statistics - Orbit types
    unclassified: Arc<AtomicUsize>,
    trapped_lost: Arc<AtomicUsize>,
    trapped_confined: Arc<AtomicUsize>,
    copassing_lost: Arc<AtomicUsize>,
    copassing_confined: Arc<AtomicUsize>,
    cupassing_lost: Arc<AtomicUsize>,
    cupassing_confined: Arc<AtomicUsize>,
    potato: Arc<AtomicUsize>,
    stagnated: Arc<AtomicUsize>,
    undefined: Arc<AtomicUsize>,
    unclassified_str: String,
    trapped_lost_str: String,
    trapped_confined_str: String,
    copassing_lost_str: String,
    copassing_confined_str: String,
    cupassing_lost_str: String,
    cupassing_confined_str: String,
    potato_str: String,
    stagnated_str: String,
    undefined_str: String,
}

impl ClassifyPbar {
    /// Initialize the pre-configured progress bar.
    pub(crate) fn new(queue: &Queue) -> Self {
        #[expect(clippy::wildcard_imports, reason = "colors")]
        use orbit_type_colors::*;
        let style = ProgressStyle::with_template(CLASSIFY_PBAR_STYLE)
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars(CLASSIFY_PROGRESS_CHARS);
        let pbar = ProgressBar::new(queue.particles.len() as u64).with_style(style);
        pbar.enable_steady_tick(Duration::from_millis(100));
        Self {
            pbar,
            length: queue.particles.len(),
            unclassified: Arc::default(),
            trapped_lost: Arc::default(),
            trapped_confined: Arc::default(),
            copassing_lost: Arc::default(),
            copassing_confined: Arc::default(),
            cupassing_lost: Arc::default(),
            cupassing_confined: Arc::default(),
            potato: Arc::default(),
            stagnated: Arc::default(),
            undefined: Arc::default(),
            unclassified_str: "Unclassified".color(UNCLASSIFIED).to_string(),
            trapped_lost_str: "Trapped-Lost".color(TRAPPED_LOST).to_string(),
            trapped_confined_str: "Trapped-Confined".color(TRAPPED_CONFINED).to_string(),
            copassing_lost_str: "CoPassing-Lost".color(COPASSING_LOST).to_string(),
            copassing_confined_str: "CoPassing-Confined".color(COPASSING_CONFINED).to_string(),
            cupassing_lost_str: "Cupassing-Lost".color(CUPASSING_LOST).to_string(),
            cupassing_confined_str: "Cupassing-Confined".color(CUPASSING_CONFINED).to_string(),
            potato_str: "Potato".color(POTATO).to_string(),
            stagnated_str: "Stagnated".color(STAGNATED).to_string(),
            undefined_str: "Undefined".color(UNDEFINED).to_string(),
        }
    }

    /// Prints an informative message before the ticking starts.
    pub(crate) fn print_prelude(&self) {
        self.pbar.println(format!(
            "🚀 Using {} threads for {} particles",
            rayon::current_num_threads(),
            self.length
        ));
        self.pbar.println("🗿 Classifying orbits");
    }

    /// Increases the wrapped pbar, as well as the live statistics.
    pub(crate) fn inc(&self, orbit_type: &OrbitType) {
        self.pbar.inc(1);
        let _: usize = match *orbit_type {
            OrbitType::Unclassified => self.unclassified.fetch_add(1, SeqCst),
            OrbitType::TrappedLost => self.trapped_lost.fetch_add(1, SeqCst),
            OrbitType::TrappedConfined => self.trapped_confined.fetch_add(1, SeqCst),
            OrbitType::CoPassingLost => self.copassing_lost.fetch_add(1, SeqCst),
            OrbitType::CoPassingConfined => self.copassing_confined.fetch_add(1, SeqCst),
            OrbitType::CuPassingLost => self.cupassing_lost.fetch_add(1, SeqCst),
            OrbitType::CuPassingConfined => self.cupassing_confined.fetch_add(1, SeqCst),
            OrbitType::Potato => self.potato.fetch_add(1, SeqCst),
            OrbitType::Stagnated => self.stagnated.fetch_add(1, SeqCst),
            OrbitType::Undefined => self.undefined.fetch_add(1, SeqCst),
            _ => unimplemented!(),
        };
    }

    /// Updates the printed live statistics.
    pub(crate) fn print_stats(&self) {
        self.pbar.set_message(format!(
            concat!(
                "======== 📊 Stats ========\n",
                "🪤 {}       = {}\n",
                "🪤 {}   = {}\n",
                "➕ {}     = {}\n",
                "➕ {} = {}\n",
                "➖ {}     = {}\n",
                "➖ {} = {}\n",
                "🥔 {}             = {}\n",
                "🛏️  {}          = {}\n",
                "❓ {}       = {}\n",
                "❔ {}          = {}\n",
            ),
            self.trapped_lost_str,
            self.trapped_lost.load(SeqCst),
            self.trapped_confined_str,
            self.trapped_confined.load(SeqCst),
            self.copassing_lost_str,
            self.copassing_lost.load(SeqCst),
            self.copassing_confined_str,
            self.copassing_confined.load(SeqCst),
            self.cupassing_lost_str,
            self.cupassing_lost.load(SeqCst),
            self.cupassing_confined_str,
            self.cupassing_confined.load(SeqCst),
            self.potato_str,
            self.potato.load(SeqCst),
            self.stagnated_str,
            self.stagnated.load(SeqCst),
            self.unclassified_str,
            self.unclassified.load(SeqCst),
            self.undefined_str,
            self.undefined.load(SeqCst),
        ));
    }

    pub(crate) fn finish(&self) {
        self.pbar.println("👌 Orbit classification done");
        self.pbar.force_draw();
        self.pbar.finish();
    }
}
