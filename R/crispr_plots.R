# crispr_plots.R — CRISPR screen plotting utilities for MOBA-seq experiments
# Reusable functions for visualizing sgRNA cell-number distributions,
# replicate correlations, ratio histograms, and low-abundance Venn diagrams.

# ---- Helper ----------------------------------------------------------------

#' Discrete colour palette
#'
#' Returns a vector of `n` colours from an RColorBrewer palette, falling back
#' to [grDevices::hcl.colors()] when RColorBrewer is not installed.
#'
#' @param n.colours Number of colours to return.
#' @param palette.name Name of an RColorBrewer palette (e.g. `"Set1"`).
#' @return A character vector of hex colour codes.
#' @export
get_palette <- function(n.colours, palette.name) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    RColorBrewer::brewer.pal(max(3, n.colours), palette.name)[seq_len(n.colours)]
  } else {
    grDevices::hcl.colors(n.colours, "Dark 3")
  }
}

# ---- plot_neg_correlation ---------------------------------------------------

#' sgRNA cell-number correlation scatter plot
#'
#' Scatter plot of one x-sample vs one or more y-samples with per-series
#' linear regression lines annotated with their equation and R-squared.
#' Works for comparing a pre-sorted (unsorted) sample against negative
#' replicates, or replicate-vs-replicate consistency checks.
#'
#' @param sample.info Path to a CSV file or a data.frame with at least
#'   columns for sgRNA id, sample id, and a numeric value.
#' @param x.sample Sample name on the x-axis.
#' @param y.samples Character vector of one or more sample names on the y-axis.
#' @param id.column Column name for the per-sgRNA identifier.
#' @param sample.column Column name holding sample names.
#' @param value.column Column name for the numeric value to correlate.
#' @param min.cells Minimum cell-number threshold applied to raw values in
#'   both axes.
#' @param min.y.cells Minimum cell-number threshold applied to y-axis values
#'   only (before log transform). Points below this in any y-sample are dropped
#'   before regression. Use to exclude low-count observations from one axis.
#' @param include.zeros If `TRUE`, guides present in `x.sample` but absent
#'   from a y-sample are kept with y set to `zero.floor.cells`.
#' @param zero.floor.cells Plotted value for absent guides (used only when
#'   `include.zeros = TRUE`).
#' @param log.transform If `TRUE`, axes are log10-transformed and regression
#'   is fit in log space.
#' @param pseudocount Value added before log10 to handle zeros.
#' @param point.size,point.alpha,line.size Aesthetic parameters.
#' @param palette.name RColorBrewer palette name.
#' @param facet.plot If `TRUE`, one panel per y-sample.
#' @param annotate.eq If `TRUE`, draw the regression equation and R-squared.
#' @param plot.title Plot title; auto-generated if `NULL`.
#' @param save.plot If `TRUE`, save to `output.file`.
#' @param output.file File path for saved plot.
#' @param plot.width,plot.height Dimensions in inches for saved plot.
#' @return A ggplot object.
#' @export
plot_neg_correlation <- function(
    sample.info,
    x.sample      = "unsorted-1",
    y.samples     = c("neg-1", "neg-2"),
    id.column     = "sgID",
    sample.column = "Sample_ID",
    value.column  = "cell_num",
    min.cells     = 0,
    min.y.cells   = 0,
    include.zeros = FALSE,
    zero.floor.cells = 0.5,
    log.transform = TRUE,
    pseudocount   = 1,
    point.size    = 0.6,
    point.alpha   = 0.25,
    line.size     = 0.9,
    palette.name  = "Set1",
    facet.plot    = FALSE,
    annotate.eq   = TRUE,
    plot.title    = NULL,
    save.plot     = FALSE,
    output.file   = "neg_correlation_scatter.pdf",
    plot.width    = 7,
    plot.height   = 6
) {
  require(tidyr)
  require(ggplot2)

  ## Read and filter
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- as.data.frame(sample.info)
  }

  keep.samples    <- c(x.sample, y.samples)
  missing.samples <- setdiff(keep.samples, unique(df[[sample.column]]))
  if (length(missing.samples) > 0) {
    stop("Sample(s) not found in '", sample.column, "': ",
         paste(missing.samples, collapse = ", "))
  }

  df <- df[df[[sample.column]] %in% keep.samples,
           c(id.column, sample.column, value.column)]

  ## Pivot wide
  wide <- tidyr::pivot_wider(
    df,
    id_cols     = dplyr::all_of(id.column),
    names_from  = dplyr::all_of(sample.column),
    values_from = dplyr::all_of(value.column)
  )
  wide <- as.data.frame(wide)

  ## Long plotting frame

  n.id <- nrow(wide)
  plot.df <- data.frame(
    sgID      = rep(wide[[id.column]], times = length(y.samples)),
    Replicate = rep(y.samples, each = n.id),
    x.value   = rep(wide[[x.sample]], times = length(y.samples)),
    y.value   = unlist(wide[, y.samples, drop = FALSE], use.names = FALSE),
    stringsAsFactors = FALSE
  )

  plot.df <- plot.df[!is.na(plot.df$x.value), ]
  if (include.zeros) {
    plot.df$y.value[is.na(plot.df$y.value)] <- zero.floor.cells
  } else {
    plot.df <- plot.df[!is.na(plot.df$y.value), ]
  }

  ## Cell-number cutoff
  if (min.cells > 0) {
    n.before <- nrow(plot.df)
    plot.df  <- plot.df[plot.df$x.value >= min.cells & plot.df$y.value >= min.cells, ]
    message(sprintf("min.cells = %g: kept %d of %d sgRNA-pairs (>= %g cells in both %s and the negative).",
                    min.cells, nrow(plot.df), n.before, min.cells, x.sample))
  }

  ## Y-axis-only cutoff (raw scale, before log)
  if (min.y.cells > 0) {
    n.before <- nrow(plot.df)
    plot.df  <- plot.df[plot.df$y.value >= min.y.cells, ]
    message(sprintf("min.y.cells = %g: kept %d of %d sgRNA-pairs (>= %g cells in y-sample).",
                    min.y.cells, nrow(plot.df), n.before, min.y.cells))
  }

  ## Log10 transform
  if (log.transform) {
    plot.df$x.value <- log10(plot.df$x.value + pseudocount)
    plot.df$y.value <- log10(plot.df$y.value + pseudocount)
    x.lab <- sprintf("log10(%s %s + %g)", x.sample, value.column, pseudocount)
    y.base <- if (length(y.samples) == 1) y.samples else "negative replicate"
    y.lab <- sprintf("log10(%s %s + %g)", y.base, value.column, pseudocount)
  } else {
    x.lab  <- paste(x.sample, value.column)
    y.base <- if (length(y.samples) == 1) y.samples else "negative replicate"
    y.lab  <- paste(y.base, value.column)
  }

  ## Per-replicate regression
  cor.stats <- do.call(rbind, lapply(y.samples, function(rep.name) {
    sub <- plot.df[plot.df$Replicate == rep.name, ]
    fit <- lm(y.value ~ x.value, data = sub)
    co  <- coef(fit)
    data.frame(
      Replicate = rep.name, n = nrow(sub),
      slope = unname(co[2]), intercept = unname(co[1]),
      r  = cor(sub$x.value, sub$y.value, method = "pearson"),
      R2 = summary(fit)$r.squared,
      stringsAsFactors = FALSE
    )
  }))

  cor.stats$eq.label  <- sprintf("y = %.3gx %+.3g,  R² = %.3f",
                                 cor.stats$slope, cor.stats$intercept, cor.stats$R2)
  cor.stats$eq.prefix <- sprintf("%s:  y = %.3gx %+.3g,  R² = %.3f",
                                 cor.stats$Replicate, cor.stats$slope,
                                 cor.stats$intercept, cor.stats$R2)
  cor.stats$Replicate <- factor(cor.stats$Replicate, levels = y.samples)

  message("Regression of ", value.column,
          if (log.transform) " (log10 scale):" else " (linear scale):")
  print(cor.stats[, c("Replicate", "n", "slope", "intercept", "r", "R2")],
        row.names = FALSE)

  ## Title
  if (is.null(plot.title)) {
    plot.title <- sprintf("sgRNA %s: %s vs %s", value.column, x.sample,
                          paste(y.samples, collapse = ", "))
  }

  ## Build plot
  plot.df$Replicate <- factor(plot.df$Replicate, levels = y.samples)
  cols <- get_palette(length(y.samples), palette.name)

  p <- ggplot(plot.df, aes(x = x.value, y = y.value, colour = Replicate)) +
    geom_point(size = point.size, alpha = point.alpha) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = line.size) +
    scale_colour_manual(values = stats::setNames(cols, y.samples)) +
    labs(x = x.lab, y = y.lab, colour = "Replicate", title = plot.title) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ## Equation annotations
  if (annotate.eq) {
    if (facet.plot) {
      p <- p +
        geom_text(data = cor.stats,
                  aes(x = -Inf, y = Inf, label = eq.label),
                  hjust = -0.05, vjust = 1.4, size = 3,
                  colour = "black", inherit.aes = FALSE)
    } else {
      x.lo <- min(plot.df$x.value); x.hi <- max(plot.df$x.value)
      y.lo <- min(plot.df$y.value); y.hi <- max(plot.df$y.value)
      step <- 0.055 * (y.hi - y.lo)
      for (i in seq_along(y.samples)) {
        lab.i <- if (length(y.samples) == 1) cor.stats$eq.label[i] else cor.stats$eq.prefix[i]
        p <- p + annotate("text",
                          x = x.lo + 0.02 * (x.hi - x.lo),
                          y = y.hi - (i - 1) * step,
                          label = lab.i, colour = cols[i],
                          hjust = 0, vjust = 1, size = 3.2)
      }
    }
  }

  ## Facet / legend
  if (facet.plot) {
    p <- p + facet_wrap(~ Replicate) + theme(legend.position = "none")
  } else if (length(y.samples) == 1) {
    p <- p + theme(legend.position = "none")
  }

  ## Save
  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  return(p)
}

# ---- plot_ratio_histogram ---------------------------------------------------

#' Ratio histogram of negative / pre-sorted cell numbers
#'
#' For each negative replicate, computes the ratio of cell number to the
#' pre-sorted (unsorted) sample for sgRNAs detected in both, and draws a
#' faceted histogram. Guides from highlight samples (e.g. positive sorts) are
#' marked as coloured points on the histogram with the lowest-ratio guides
#' listed in a corner annotation box.
#'
#' @inheritParams plot_neg_correlation
#' @param presorted.sample Name of the denominator (pre-sorted / unsorted)
#'   sample.
#' @param neg.samples Character vector of negative-sort sample names (one
#'   histogram panel each).
#' @param highlight.samples Character vector of samples whose guides are
#'   marked on the histogram.
#' @param log.ratio If `TRUE`, x-axis is log10(ratio).
#' @param bins Number of histogram bins.
#' @param n.label Number of lowest-x highlight guides listed per panel.
#' @param label.digits Decimal places shown for the listed values.
#' @param point.size,label.size,line.height Aesthetic parameters.
#' @param hist.fill,hist.colour Histogram bar colours.
#' @return A ggplot object.
#' @export
plot_ratio_histogram <- function(
    sample.info,
    presorted.sample  = "unsorted-1",
    neg.samples       = c("neg-1", "neg-2"),
    highlight.samples = c("pos-1", "pos-2"),
    id.column     = "sgID",
    sample.column = "Sample_ID",
    value.column  = "cell_num",
    log.ratio     = TRUE,
    bins          = 50,
    n.label       = 10,
    label.digits  = 2,
    point.size    = 1.8,
    label.size    = 2.6,
    line.height   = 1.25,
    palette.name  = "Set1",
    hist.fill     = "grey80",
    hist.colour   = "grey40",
    plot.title    = NULL,
    save.plot     = FALSE,
    output.file   = "ratio_histogram.pdf",
    plot.width    = 10,
    plot.height   = 7
) {
  require(ggplot2)

  ## Read and validate
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- as.data.frame(sample.info)
  }
  needed <- c(presorted.sample, neg.samples, highlight.samples)
  missing.samples <- setdiff(needed, unique(df[[sample.column]]))
  if (length(missing.samples) > 0) {
    stop("Sample(s) not found in '", sample.column, "': ",
         paste(missing.samples, collapse = ", "))
  }
  id <- id.column; sc <- sample.column; vc <- value.column

  ## Pre-sorted denominator
  ns <- df[df[[sc]] == presorted.sample & df[[vc]] > 0, c(id, vc)]
  names(ns) <- c("sgID", "presorted")

  ## Ratio per negative replicate
  ratio.df <- do.call(rbind, lapply(neg.samples, function(neg) {
    sub <- df[df[[sc]] == neg & df[[vc]] > 0, c(id, vc)]
    names(sub) <- c("sgID", "neg")
    m <- merge(sub, ns, by = "sgID")
    if (nrow(m) == 0) return(NULL)
    m$Replicate <- neg
    m$ratio     <- m$neg / m$presorted
    m
  }))
  ratio.df$Replicate <- factor(ratio.df$Replicate, levels = neg.samples)
  ratio.df$x <- if (log.ratio) log10(ratio.df$ratio) else ratio.df$ratio

  ## Highlight guides
  hl <- df[df[[sc]] %in% highlight.samples & df[[vc]] > 0, c(id, sc, vc)]
  names(hl) <- c("sgID", "origin", "value")
  hl <- hl[order(hl$sgID, -hl$value), ]
  hl <- hl[!duplicated(hl$sgID), c("sgID", "origin")]
  hl$origin <- factor(hl$origin, levels = highlight.samples)

  hl.pos <- merge(ratio.df[, c("sgID", "Replicate", "x")], hl, by = "sgID")

  ## Annotation box
  box.df <- do.call(rbind, lapply(split(hl.pos, hl.pos$Replicate), function(d) {
    if (nrow(d) == 0) return(NULL)
    d <- d[order(d$x), ]
    d <- d[seq_len(min(nrow(d), n.label)), ]
    d$row <- seq_len(nrow(d))
    d
  }))
  box.df$box.label <- sprintf(paste0("%s  %.", label.digits, "f"),
                              box.df$sgID, box.df$x)
  box.df$vpos <- box.df$row * line.height

  ## Report
  message(sprintf("Ratio = negative / %s%s", presorted.sample,
                  if (log.ratio) " (log10)" else ""))
  rep.summary <- as.data.frame(table(ratio.df$Replicate))
  names(rep.summary) <- c("Replicate", "n_sgRNA")
  hl.summary <- as.data.frame(table(hl.pos$Replicate))
  names(hl.summary) <- c("Replicate", "n_highlight_points")
  print(merge(rep.summary, hl.summary, by = "Replicate"), row.names = FALSE)

  ## Title / axis
  if (is.null(plot.title)) {
    plot.title <- sprintf("Negative / %s cell-number ratio per sgRNA\nmarked: %s",
                          presorted.sample, paste(highlight.samples, collapse = ", "))
  }
  x.lab <- if (log.ratio) {
    sprintf("log10(negative / %s cell number)", presorted.sample)
  } else {
    sprintf("negative / %s cell number", presorted.sample)
  }

  ## Build plot
  p <- ggplot(ratio.df, aes(x = x)) +
    geom_histogram(bins = bins, fill = hist.fill, colour = hist.colour, linewidth = 0.2) +
    geom_point(data = hl.pos, aes(x = x, y = 0, colour = origin),
               size = point.size, inherit.aes = FALSE) +
    geom_text(data = box.df,
              aes(x = Inf, y = Inf, label = box.label, colour = origin),
              hjust = 1.02, vjust = box.df$vpos, size = label.size,
              inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(~ Replicate, scales = "free_y") +
    scale_colour_brewer(palette = palette.name, drop = FALSE) +
    labs(x = x.lab, y = "Number of sgRNAs",
         colour = "Sample of origin", title = plot.title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  ## Save
  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  return(p)
}

# ---- plot_neg_hits_scatter --------------------------------------------------

#' Faceted scatter with positive hits highlighted
#'
#' Faceted scatter of a pre-sorted sample (x) vs each negative replicate (y),
#' with all points in grey and each panel's corresponding positive-sample hits
#' highlighted in colour and labelled. Correspondence defaults to suffix
#' matching (`sub("neg", "pos", ...)`). A horizontal cutoff line and black
#' rings around below-cutoff hits are drawn by default.
#'
#' @inheritParams plot_neg_correlation
#' @param neg.samples Character vector of negative sample names.
#' @param hit.map Named character vector mapping each negative sample to its
#'   corresponding positive sample.
#' @param spike.sample Optional spike-in sample name to overlay; `NULL` to
#'   disable.
#' @param cutoff.cells Horizontal cutoff line position in cells.
#' @param draw.cutoff Whether to draw the cutoff line.
#' @param border.below Whether to add a black ring around hits below cutoff.
#' @param label.hits Which hits to label: `"below"`, `"all"`, or `"none"`.
#' @param label.markers Extra sgIDs to label if present as a hit.
#' @param base.colour Colour for non-hit background points.
#' @param hit.size,label.size,spike.colour,spike.shape,spike.size Aesthetics.
#' @return A ggplot object.
#' @export
plot_neg_hits_scatter <- function(
    sample.info,
    x.sample      = "unsorted-1",
    neg.samples   = c("neg-1", "neg-2"),
    hit.map       = NULL,
    spike.sample  = NULL,
    id.column     = "sgID",
    sample.column = "Sample_ID",
    value.column  = "cell_num",
    log.transform = TRUE,
    pseudocount   = 0,
    min.cells     = 0,
    include.zeros = FALSE,
    zero.floor.cells = 0.5,
    cutoff.cells  = 2,
    draw.cutoff   = TRUE,
    border.below  = TRUE,
    label.hits    = "below",
    label.markers = NULL,
    base.colour   = "grey75",
    palette.name  = "Set1",
    point.size    = 0.6,
    point.alpha   = 0.30,
    hit.size      = 1.9,
    label.size    = 2.6,
    spike.colour  = "gray35",
    spike.shape   = 17,
    spike.size    = 2.4,
    plot.title    = NULL,
    save.plot     = FALSE,
    output.file   = "neg_hits_scatter.pdf",
    plot.width    = 9,
    plot.height   = 7
) {
  require(tidyr)
  require(ggplot2)

  ## Resolve neg -> pos correspondence
  if (is.null(hit.map)) hit.map <- stats::setNames(sub("^neg", "pos", neg.samples), neg.samples)
  pos.samples <- unname(hit.map)

  ## Read and validate
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else df <- as.data.frame(sample.info)
  needed <- c(x.sample, neg.samples, pos.samples, spike.sample)
  missing.samples <- setdiff(needed, unique(df[[sample.column]]))
  if (length(missing.samples) > 0)
    stop("Sample(s) not found in '", sample.column, "': ", paste(missing.samples, collapse = ", "))
  id <- id.column; sc <- sample.column; vc <- value.column

  ## Wide
  sub.df <- df[df[[sc]] %in% c(x.sample, neg.samples), c(id, sc, vc)]
  wide <- as.data.frame(tidyr::pivot_wider(
    sub.df, id_cols = dplyr::all_of(id),
    names_from = dplyr::all_of(sc), values_from = dplyr::all_of(vc)))

  ## Long plotting frame
  n.id <- nrow(wide)
  plot.df <- data.frame(
    sgID      = rep(wide[[id]], times = length(neg.samples)),
    Replicate = rep(neg.samples, each = n.id),
    x.value   = rep(wide[[x.sample]], times = length(neg.samples)),
    y.value   = unlist(wide[, neg.samples, drop = FALSE], use.names = FALSE),
    stringsAsFactors = FALSE)

  plot.df <- plot.df[!is.na(plot.df$x.value), ]
  plot.df$is.zero <- is.na(plot.df$y.value)

  if (include.zeros) {
    plot.df$y.value[plot.df$is.zero] <- zero.floor.cells
  } else {
    plot.df <- plot.df[!plot.df$is.zero, ]
  }

  if (min.cells > 0) {
    keep <- plot.df$x.value >= min.cells & (plot.df$is.zero | plot.df$y.value >= min.cells)
    plot.df <- plot.df[keep, ]
  }

  ## Flag corresponding positive hits
  plot.df$Hit.sample <- NA_character_
  for (neg in neg.samples) {
    pos     <- hit.map[[neg]]
    hit.ids <- unique(df[[id]][df[[sc]] == pos & df[[vc]] > 0])
    plot.df$Hit.sample[plot.df$Replicate == neg & plot.df$sgID %in% hit.ids] <- pos
  }
  plot.df$Replicate  <- factor(plot.df$Replicate,  levels = neg.samples)
  plot.df$Hit.sample <- factor(plot.df$Hit.sample, levels = pos.samples)

  ## Transform
  cutoff.y <- if (log.transform) log10(cutoff.cells + pseudocount) else cutoff.cells
  if (log.transform) {
    plot.df$x.value <- log10(plot.df$x.value + pseudocount)
    plot.df$y.value <- log10(plot.df$y.value + pseudocount)
    x.lab <- sprintf("log10(%s %s)", x.sample, value.column)
    y.lab <- sprintf("log10(negative replicate %s)", value.column)
  } else {
    x.lab <- paste(x.sample, value.column); y.lab <- paste("negative replicate", value.column)
  }

  hits.df  <- plot.df[!is.na(plot.df$Hit.sample), ]

  ## Spike-ins
  if (!is.null(spike.sample)) {
    spike.ids <- unique(df[[id]][df[[sc]] == spike.sample & df[[vc]] > 0])
    spike.df  <- plot.df[plot.df$sgID %in% spike.ids, ]
  } else spike.df <- plot.df[0, ]

  ## Hits to label / ring
  border.df <- if (border.below) hits.df[hits.df$y.value < cutoff.y, ] else hits.df[0, ]
  label.df  <- switch(label.hits,
                      below = hits.df[hits.df$y.value < cutoff.y, ],
                      all   = hits.df,
                      none  = hits.df[0, ],
                      hits.df[hits.df$y.value < cutoff.y, ])
  marker.df <- hits.df[!(hits.df$sgID %in% label.df$sgID) & hits.df$sgID %in% label.markers, ]

  ## Report
  rep.tab <- as.data.frame(table(Replicate = hits.df$Replicate))
  names(rep.tab)[2] <- "n_hits"
  bel.tab <- as.data.frame(table(Replicate = border.df$Replicate)); names(bel.tab)[2] <- "n_below_cutoff"
  print(merge(rep.tab, bel.tab, by = "Replicate"), row.names = FALSE)

  ## Title
  if (is.null(plot.title))
    plot.title <- sprintf("%s vs negative replicates (cutoff = %g cells)", x.sample, cutoff.cells)

  ## Build plot
  p <- ggplot(plot.df, aes(x = x.value, y = y.value)) +
    geom_point(colour = base.colour, size = point.size, alpha = point.alpha) +
    geom_point(data = hits.df, aes(colour = Hit.sample), size = hit.size) +
    facet_wrap(~ Replicate) +
    scale_colour_brewer(palette = palette.name, drop = FALSE) +
    labs(x = x.lab, y = y.lab, colour = "Positive sample (hits)", title = plot.title,
         caption = sprintf("Dashed line: %g cells; black ring: hits below cutoff%s",
                           cutoff.cells,
                           if (include.zeros) "; zeros placed at floor" else "")) +
    theme_bw(base_size = 11) + theme(legend.position = "right")

  if (nrow(spike.df) > 0)
    p <- p + geom_point(data = spike.df, aes(x = x.value, y = y.value),
                        colour = spike.colour, shape = spike.shape, size = spike.size,
                        inherit.aes = FALSE)

  if (draw.cutoff)
    p <- p + geom_hline(yintercept = cutoff.y, linetype = "dashed", colour = "grey30", linewidth = 0.5)

  if (nrow(border.df))
    p <- p + geom_point(data = border.df, aes(x = x.value, y = y.value),
                        shape = 1, colour = "black", size = hit.size + 1.3,
                        stroke = 0.6, inherit.aes = FALSE)

  add_lab <- function(plot, d) {
    if (nrow(d) == 0) return(plot)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot + ggrepel::geom_text_repel(
        data = d, aes(label = sgID, colour = Hit.sample),
        size = label.size, show.legend = FALSE, max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.2, box.padding = 0.3)
    } else {
      plot + geom_text(data = d, aes(label = sgID, colour = Hit.sample),
                       size = label.size, vjust = -0.6,
                       show.legend = FALSE, check_overlap = TRUE)
    }
  }
  p <- add_lab(p, label.df)
  p <- add_lab(p, marker.df)

  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  return(p)
}

# ---- plot_low_abundance_venn ------------------------------------------------

#' Venn diagram of low-abundance sgRNAs across negative samples
#'
#' Draws a Venn diagram of sgRNAs (or genes) whose cell number falls below a
#' threshold in each negative sample. Optionally overlays positive-hit
#' membership and prints a table of shared hits.
#'
#' @inheritParams plot_neg_correlation
#' @param neg.samples Character vector of negative sample names.
#' @param hit.samples Optional character vector of positive samples whose
#'   guides are overlaid on the Venn.
#' @param log.threshold sgRNAs with log10(cell_num) below this are counted as
#'   low-abundance.
#' @param include.zeros If `TRUE`, sgRNAs absent from a negative but present
#'   in `reference.sample` are counted as low-abundance (0 cells).
#' @param reference.sample Universe sample for the zero rule.
#' @param gene.level If `TRUE`, collapse sgRNAs to genes using `sg.pattern`.
#' @param sg.pattern Regex to strip the sgRNA suffix for gene-level analysis.
#' @param annotate `"count"` to show hit counts per Venn region, or `"none"`.
#' @param name.shared.min If set, table hits low in at least this many
#'   negatives.
#' @param name.max Maximum rows in the corner table.
#' @return A list with `plot`, `sets`, and `hits` (invisibly).
#' @export
plot_low_abundance_venn <- function(
    sample.info,
    neg.samples   = c("neg-1", "neg-2"),
    hit.samples   = NULL,
    id.column     = "sgID",
    sample.column = "Sample_ID",
    value.column  = "cell_num",
    log.threshold = log10(2),
    include.zeros = TRUE,
    reference.sample = "unsorted-1",
    gene.level    = FALSE,
    sg.pattern    = "_sg[0-9]+$",
    annotate      = "count",
    name.shared.min = NULL,
    name.max      = 12,
    palette.name  = "Set1",
    plot.title    = NULL,
    save.plot     = FALSE,
    output.file   = "low_abundance_venn.pdf",
    plot.width    = 7,
    plot.height   = 6
) {
  require(ggplot2)

  ## Read and validate
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else df <- as.data.frame(sample.info)
  needed <- c(neg.samples, hit.samples, if (include.zeros) reference.sample)
  missing.samples <- setdiff(needed, unique(df[[sample.column]]))
  if (length(missing.samples) > 0)
    stop("Sample(s) not found in '", sample.column, "': ", paste(missing.samples, collapse = ", "))
  sc <- sample.column; vc <- value.column; idc <- id.column
  to.feat <- function(x) if (gene.level) sub(sg.pattern, "", x) else x

  ## Low-abundance feature set per negative
  universe <- if (include.zeros) unique(to.feat(df[[idc]][df[[sc]] == reference.sample & df[[vc]] > 0])) else NULL
  set.list <- lapply(neg.samples, function(s) {
    sub <- df[df[[sc]] == s & df[[vc]] > 0, ]
    low.present <- unique(to.feat(sub[[idc]][log10(sub[[vc]]) < log.threshold]))
    if (include.zeros) {
      present.here <- unique(to.feat(sub[[idc]]))
      union(low.present, setdiff(universe, present.here))
    } else low.present
  })
  names(set.list) <- neg.samples

  message(sprintf("Low-abundance %s with log10(%s) < %g:",
                  if (gene.level) "genes" else "sgRNAs", vc, log.threshold))
  for (s in neg.samples) message(sprintf("  %-9s %d", s, length(set.list[[s]])))
  message(sprintf("  union %d | shared %d",
                  length(Reduce(union, set.list)), length(Reduce(intersect, set.list))))

  ## Map positive hits -> region
  hit.table <- NULL
  if (!is.null(hit.samples)) {
    hit.table <- do.call(rbind, lapply(hit.samples, function(p) {
      genes <- unique(to.feat(df[[idc]][df[[sc]] == p & df[[vc]] > 0]))
      if (length(genes) == 0) return(NULL)
      region <- vapply(genes, function(g) {
        mem <- neg.samples[vapply(set.list, function(s) g %in% s, logical(1))]
        if (length(mem) == 0) "none" else paste(mem, collapse = "/")
      }, character(1))
      data.frame(gene = genes, sample = p, region = region,
                 n.neg = lengths(lapply(genes, function(g)
                   neg.samples[vapply(set.list, function(s) g %in% s, logical(1))])),
                 stringsAsFactors = FALSE)
    }))
    hit.table$sample <- factor(hit.table$sample, levels = hit.samples)
    cat("\nPositive hits by negative-low-abundance region:\n")
    print(hit.table[order(hit.table$sample, -hit.table$n.neg, hit.table$gene),
                    c("sample", "gene", "n.neg", "region")], row.names = FALSE)
    n.none <- sum(hit.table$region == "none")
    if (n.none > 0)
      message(sprintf("\n%d hit feature(s) are NOT low in any negative (cannot sit in a Venn region).",
                      n.none))
    if (!is.null(name.shared.min)) {
      sh <- hit.table[hit.table$n.neg >= name.shared.min, ]
      cat(sprintf("\nPositive hits low in >= %d negatives (%d features):\n",
                  name.shared.min, nrow(sh)))
      print(sh[order(sh$sample, -sh$n.neg, sh$gene), c("sample", "gene", "n.neg", "region")],
            row.names = FALSE)
    }
  }

  ## Title
  if (is.null(plot.title))
    plot.title <- sprintf("%s with log10(%s) < %g across negative samples",
                          if (gene.level) "Genes" else "sgRNAs", vc, log.threshold)

  ## Draw Venn
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    warning("Install 'ggVennDiagram' to draw the Venn. Returning sets and hit table only.")
    return(invisible(list(plot = NULL, sets = set.list, hits = hit.table)))
  }

  p <- ggVennDiagram::ggVennDiagram(set.list, label = "count", label_alpha = 0) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF", guide = "none") +
    ggplot2::labs(title = plot.title)

  if (!is.null(hit.samples)) {
    p <- p + ggplot2::scale_colour_brewer(palette = palette.name, name = "Positive sample") +
      ggplot2::theme(legend.position = "right")

    if (annotate == "count" && requireNamespace("sf", quietly = TRUE)) {
      tryCatch({
        pdata <- ggVennDiagram::process_data(ggVennDiagram::Venn(set.list))
        rl    <- ggVennDiagram::venn_regionlabel(pdata)
        xy    <- sf::st_coordinates(sf::st_geometry(rl))
        key.col <- intersect(c("name", "id"), names(rl))[1]
        rl$key  <- as.character(rl[[key.col]])
        rl$cx <- xy[, 1]; rl$cy <- xy[, 2]
        norm.key <- function(k) {
          parts <- unlist(strsplit(k, "/", fixed = TRUE))
          paste(neg.samples[neg.samples %in% parts], collapse = "/")
        }
        rl$norm <- vapply(rl$key, norm.key, character(1))

        cnt <- as.data.frame(table(region = hit.table$region, sample = hit.table$sample))
        cnt <- cnt[cnt$Freq > 0 & cnt$region != "none", ]
        cnt <- merge(cnt, rl[, c("norm", "cx", "cy")], by.x = "region", by.y = "norm")
        cnt <- cnt[order(cnt$region, cnt$sample), ]
        cnt$k  <- ave(seq_len(nrow(cnt)), cnt$region, FUN = seq_along)
        cnt$ny <- cnt$cy - (cnt$k - 1) * 0.022 + 0.011 * (ave(cnt$k, cnt$region, FUN = max) - 1)
        cnt$lab <- sprintf("%s: %d", cnt$sample, cnt$Freq)
        p <<- p + ggplot2::geom_text(
          data = cnt, ggplot2::aes(x = cx, y = ny, label = lab, colour = sample),
          size = 2.8, fontface = "bold", show.legend = TRUE, inherit.aes = FALSE)
      }, error = function(e)
        message("On-Venn count placement failed (", conditionMessage(e),
                "); region counts omitted, hit table still printed."))
    }

    if (!is.null(name.shared.min)) {
      sh <- hit.table[hit.table$n.neg >= name.shared.min, ]
      if (nrow(sh) > 0) {
        sh <- sh[order(sh$sample, -sh$n.neg, sh$gene), ]
        n.total <- nrow(sh)
        if (n.total > name.max) sh <- sh[seq_len(name.max), ]
        sh$line <- sprintf("%s  (%d)", sh$gene, sh$n.neg)
        rows  <- c(sprintf("low in ≥ %d negatives:", name.shared.min), sh$line,
                   if (n.total > name.max) sprintf("+%d more (see $hits)", n.total - name.max))
        cols  <- c("grey15", as.character(sh$sample),
                   if (n.total > name.max) "grey40")
        n.row <- length(rows)
        tab.df <- data.frame(
          x = Inf, y = -Inf, label = rows, col = cols,
          vj = -(n.row - seq_len(n.row)) * 1.05,
          face = c("bold", rep("plain", n.row - 1)),
          stringsAsFactors = FALSE)
        pal <- if (requireNamespace("RColorBrewer", quietly = TRUE))
          stats::setNames(RColorBrewer::brewer.pal(max(3, length(hit.samples)), palette.name)[seq_along(hit.samples)],
                          hit.samples) else NULL
        tab.df$col <- ifelse(tab.df$col %in% names(pal), pal[tab.df$col], tab.df$col)
        p <- p + ggplot2::geom_text(
          data = tab.df, ggplot2::aes(x = x, y = y, label = label),
          colour = tab.df$col, fontface = tab.df$face, vjust = tab.df$vj,
          hjust = 1.02, size = 2.5, inherit.aes = FALSE)
      }
    }
  }

  print(p)
  if (save.plot) {
    ggplot2::ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }
  invisible(list(plot = p, sets = set.list, hits = hit.table))
}

# ---- plot_ratio_histogram_hits ----------------------------------------------

#' Ratio histogram with corresponding positive hits and statistics
#'
#' Per-replicate histogram of (negative / pre-sorted) cell-number ratio, where
#' each panel shows only its corresponding positive sample's hits. For each
#' panel it finds the lowest-x hit, z-scores it against the ratio distribution,
#' and reports the fraction of negative cells at or below that threshold.
#'
#' @inheritParams plot_neg_correlation
#' @param presorted.sample Name of the denominator sample.
#' @param neg.samples Character vector of negative sample names.
#' @param hit.map Named character vector mapping neg -> pos samples.
#' @param spike.sample Optional spike-in sample; `NULL` to disable.
#' @param show.stats If `TRUE`, annotate per-panel statistics.
#' @param log.ratio If `TRUE`, x-axis is log10(ratio).
#' @param z.tail `"lower"` for one-sided or `"two"` for two-sided p-value.
#' @param bold.quantile Bottom quantile of the ratio distribution to bold.
#' @param label.markers sgIDs to label in a distinct colour.
#' @param marker.colour Colour for marker labels.
#' @param bins Number of histogram bins.
#' @param label.mode `"list"`, `"points"`, or `"both"`.
#' @param include.zeros If `TRUE`, show guides absent from negative as zeros.
#' @param zero.x X-position for the zero bar.
#' @param zero.bar If `TRUE`, draw a capped bar for zero-count guides.
#' @return A ggplot object.
#' @export
plot_ratio_histogram_hits <- function(
    sample.info,
    presorted.sample = "unsorted-1",
    neg.samples      = c("neg-1", "neg-2"),
    hit.map          = NULL,
    spike.sample     = NULL,
    id.column     = "sgID",
    sample.column = "Sample_ID",
    show.stats    = TRUE,
    value.column  = "cell_num",
    log.ratio     = TRUE,
    z.tail        = "lower",
    bold.quantile = 0.05,
    label.markers = NULL,
    marker.colour = "black",
    bins          = 50,
    list.size     = 2.6,
    stat.size     = 2.9,
    point.size    = 1.9,
    spike.colour  = "gray35",
    spike.shape   = 17,
    spike.size    = 2.3,
    label.mode    = "both",
    palette.name  = "Set1",
    hist.fill     = "grey80",
    hist.colour   = "grey45",
    plot.title    = NULL,
    save.plot     = FALSE,
    output.file   = "ratio_histogram_hits.pdf",
    plot.width    = 11,
    plot.height   = 7,
    include.zeros = TRUE,
    zero.x        = -3,
    zero.bar      = TRUE
) {
  require(ggplot2)

  ## Resolve neg -> pos correspondence
  if (is.null(hit.map)) hit.map <- stats::setNames(sub("^neg", "pos", neg.samples), neg.samples)
  pos.samples <- unname(hit.map)

  ## Read and validate
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else df <- as.data.frame(sample.info)
  needed <- c(presorted.sample, neg.samples, pos.samples, spike.sample)
  missing.samples <- setdiff(needed, unique(df[[sample.column]]))
  if (length(missing.samples) > 0)
    stop("Sample(s) not found in '", sample.column, "': ", paste(missing.samples, collapse = ", "))
  id <- id.column; sc <- sample.column; vc <- value.column

  ## Pre-sorted denominator
  ns <- df[df[[sc]] == presorted.sample & df[[vc]] > 0, c(id, vc)]
  names(ns) <- c("sgID", "presorted")

  ## Ratio per negative replicate
  ratio.df <- do.call(rbind, lapply(neg.samples, function(neg) {
    sub <- df[df[[sc]] == neg & df[[vc]] > 0, c(id, vc)]
    names(sub) <- c("sgID", "neg")
    m <- merge(sub, ns, by = "sgID")
    if (nrow(m) == 0) return(NULL)
    m$Replicate <- neg; m$ratio <- m$neg / m$presorted; m
  }))
  ratio.df$x <- if (log.ratio) log10(ratio.df$ratio) else ratio.df$ratio
  ratio.df$Replicate <- factor(ratio.df$Replicate, levels = neg.samples)

  group.levels <- c(pos.samples, if (!is.null(spike.sample)) spike.sample)

  ## Corresponding positive hits
  hits.df <- do.call(rbind, lapply(neg.samples, function(neg) {
    hit.ids <- unique(df[[id]][df[[sc]] == hit.map[[neg]] & df[[vc]] > 0])
    d <- ratio.df[ratio.df$Replicate == neg & ratio.df$sgID %in% hit.ids, ]
    if (nrow(d) == 0) return(NULL)
    d$Group <- hit.map[[neg]]; d
  }))
  hits.df$Group <- factor(hits.df$Group, levels = group.levels)

  ## Zeros
  zero.df <- ratio.df[0, ]
  zero.hits.df <- hits.df[0, ]
  if (include.zeros) {
    pres.ids <- ns$sgID
    zero.df <- do.call(rbind, lapply(neg.samples, function(neg) {
      present.neg <- unique(df[[id]][df[[sc]] == neg & df[[vc]] > 0])
      z.ids <- setdiff(pres.ids, present.neg)
      if (length(z.ids) == 0) return(NULL)
      data.frame(sgID = z.ids, neg = 0, presorted = ns$presorted[match(z.ids, pres.ids)],
                 Replicate = neg, ratio = 0, x = zero.x, stringsAsFactors = FALSE)
    }))
    if (is.null(zero.df)) zero.df <- ratio.df[0, ]
    else zero.df$Replicate <- factor(zero.df$Replicate, levels = neg.samples)
    zero.hits.df <- do.call(rbind, lapply(neg.samples, function(neg) {
      hit.ids <- unique(df[[id]][df[[sc]] == hit.map[[neg]] & df[[vc]] > 0])
      z <- zero.df[zero.df$Replicate == neg & zero.df$sgID %in% hit.ids, ]
      if (nrow(z) == 0) return(NULL)
      z$Group <- hit.map[[neg]]; z
    }))
    if (is.null(zero.hits.df)) zero.hits.df <- hits.df[0, ]
    else zero.hits.df$Group <- factor(zero.hits.df$Group, levels = group.levels)
  }

  hist.data <- ratio.df
  hits.show <- hits.df
  if (include.zeros && nrow(zero.hits.df)) hits.show <- rbind(hits.df, zero.hits.df)

  ## Spike-ins
  if (!is.null(spike.sample)) {
    spike.ids <- unique(df[[id]][df[[sc]] == spike.sample & df[[vc]] > 0])
    spike.df  <- ratio.df[ratio.df$sgID %in% spike.ids, ]
    if (nrow(spike.df) > 0) spike.df$Group <- factor(spike.sample, levels = group.levels)
  } else {
    spike.df <- ratio.df[0, ]; spike.df$Group <- factor(character(0), levels = group.levels)
  }

  ## Bar-top y for each hit
  hbar <- ggplot2::ggplot_build(
    ggplot(hist.data, aes(x = x)) + geom_histogram(bins = bins) +
      facet_wrap(~ Replicate, scales = "free_y"))$data[[1]]
  hbar$Replicate <- neg.samples[as.integer(hbar$PANEL)]
  max.count <- tapply(hbar$count, hbar$Replicate, max)
  bw        <- stats::median(hbar$xmax - hbar$xmin)
  bar.top <- function(xv, repv) {
    b <- hbar[hbar$Replicate == repv, ]
    i <- which(xv >= b$xmin & xv <= b$xmax)
    if (length(i) == 0) NA_real_ else b$count[i[1]]
  }
  if (nrow(hits.show)) {
    hits.show$ytop <- mapply(bar.top, hits.show$x, as.character(hits.show$Replicate))
    if (include.zeros) {
      zr <- hits.show$neg == 0
      hits.show$ytop[zr] <- if (zero.bar) as.numeric(max.count[as.character(hits.show$Replicate[zr])]) else 0
    }
  }

  ## Capped zero bar
  zerobar.df <- NULL
  if (include.zeros && zero.bar && nrow(zero.df)) {
    zn <- table(factor(zero.df$Replicate, levels = neg.samples))
    zerobar.df <- data.frame(
      Replicate = factor(neg.samples, levels = neg.samples),
      x = zero.x, height = as.numeric(max.count[neg.samples]),
      n = as.integer(zn[neg.samples]), stringsAsFactors = FALSE)
    zerobar.df <- zerobar.df[!is.na(zerobar.df$height) & !is.na(zerobar.df$n) & zerobar.df$n > 0, ]
    zerobar.df$nlab <- paste0("n=", format(zerobar.df$n, big.mark = ",", trim = TRUE))
  }

  ## Per-panel statistics
  stat.df <- do.call(rbind, lapply(neg.samples, function(neg) {
    panel <- ratio.df[ratio.df$Replicate == neg, ]
    ph    <- hits.df[hits.df$Replicate == neg, ]
    mu <- mean(panel$x); sigma <- sd(panel$x); total <- sum(panel$neg)
    if (nrow(ph) == 0)
      return(data.frame(Replicate = neg, low.sgID = NA, threshold = NA, z = NA, p = NA,
                        below.cells = NA, total.cells = total, frac.below = NA,
                        stringsAsFactors = FALSE))
    i.min <- which.min(ph$x); thr <- ph$x[i.min]
    z <- (thr - mu) / sigma
    p <- if (z.tail == "two") 2 * pnorm(-abs(z)) else pnorm(z)
    below <- sum(panel$neg[panel$x <= thr])
    data.frame(Replicate = neg, low.sgID = ph$sgID[i.min], threshold = thr, z = z, p = p,
               below.cells = below, total.cells = total, frac.below = below / total,
               stringsAsFactors = FALSE)
  }))
  stat.df$Replicate <- factor(stat.df$Replicate, levels = neg.samples)
  comma <- function(v) format(round(v), big.mark = ",", trim = TRUE)
  stat.df$label <- sprintf(
    "lowest %s: %s\nx = %.2f,  z = %.2f,  p = %.2g\ncells ≤ x: %s / %s  (%.1f%%)",
    hit.map[as.character(stat.df$Replicate)], stat.df$low.sgID, stat.df$threshold,
    stat.df$z, stat.df$p, comma(stat.df$below.cells), comma(stat.df$total.cells),
    100 * stat.df$frac.below)

  ## Bottom-quantile cutoff
  qthr.map <- tapply(ratio.df$x, ratio.df$Replicate,
                     function(v) stats::quantile(v, bold.quantile, names = FALSE))
  stat.df$q.thr <- qthr.map[as.character(stat.df$Replicate)]

  ## Corner list
  header.df <- data.frame(
    Replicate = factor(neg.samples, levels = neg.samples),
    htext = paste0(hit.map[neg.samples], " hits (low→high)"), stringsAsFactors = FALSE)
  entry.df <- do.call(rbind, lapply(neg.samples, function(neg) {
    d <- hits.show[hits.show$Replicate == neg, ]
    if (nrow(d) == 0) return(NULL)
    d <- d[order(d$x), ]
    d$is.bold   <- d$x <= qthr.map[[neg]]
    d$is.marker <- d$sgID %in% label.markers
    d$row      <- seq_len(nrow(d))
    d$entry    <- ifelse(d$neg == 0, sprintf("%s  (0)", d$sgID),
                         sprintf("%s  %.2f", d$sgID, d$x))
    d$fontface <- ifelse(d$is.bold, "bold", "plain")
    d$col      <- ifelse(d$is.marker, marker.colour, "grey25")
    d
  }))
  entry.df$Replicate <- factor(entry.df$Replicate, levels = neg.samples)
  entry.df$vpos <- 2.4 + entry.df$row * 1.15

  ## Report
  if (show.stats)
    print(stat.df[, c("Replicate", "low.sgID", "threshold", "z", "p",
                      "below.cells", "total.cells", "frac.below")], row.names = FALSE)

  ## Colour / shape scales
  pos.cols <- if (requireNamespace("RColorBrewer", quietly = TRUE))
    RColorBrewer::brewer.pal(max(3, length(pos.samples)), palette.name)[seq_along(pos.samples)]
  else grDevices::hcl.colors(length(pos.samples), "Dark 3")
  grp.cols   <- stats::setNames(c(pos.cols, if (!is.null(spike.sample)) spike.colour), group.levels)
  grp.shapes <- stats::setNames(c(rep(16, length(pos.samples)),
                                   if (!is.null(spike.sample)) spike.shape), group.levels)

  ## Title / axis
  if (is.null(plot.title))
    plot.title <- sprintf("Negative / %s cell-number ratio - corresponding positive hits",
                          presorted.sample)
  x.lab <- if (log.ratio) sprintf("log10(negative / %s cell number)", presorted.sample)
  else sprintf("negative / %s cell number", presorted.sample)

  ## Build plot
  p <- ggplot(hist.data, aes(x = x)) +
    geom_histogram(bins = bins, fill = hist.fill, colour = hist.colour, linewidth = 0.2) +
    geom_vline(data = stat.df, aes(xintercept = q.thr),
               linetype = "dotted", colour = "steelblue", linewidth = 0.4) +
    geom_point(data = hits.show, aes(x = x, y = ytop, colour = Group, shape = Group),
               size = point.size, inherit.aes = FALSE) +
    facet_wrap(~ Replicate, scales = "free_y") +
    scale_colour_manual(values = grp.cols, drop = FALSE) +
    scale_shape_manual(values = grp.shapes, drop = FALSE) +
    labs(x = x.lab, y = "Number of sgRNAs", colour = "Sample", shape = "Sample",
         title = plot.title,
         caption = sprintf("Dotted: bottom %.0f%% of ratio (bold)%s", 100 * bold.quantile,
                           if (show.stats) "; dashed: lowest corresponding hit" else "")) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme_bw(base_size = 11) + theme(legend.position = "right")

  if (nrow(spike.df) > 0)
    p <- p + geom_point(data = spike.df, aes(x = x, y = 0, colour = Group, shape = Group),
                        size = spike.size, inherit.aes = FALSE)

  if (!is.null(zerobar.df)) {
    p <- p +
      geom_col(data = zerobar.df, aes(x = x, y = height), width = bw,
               fill = "grey55", colour = hist.colour, linewidth = 0.2, inherit.aes = FALSE) +
      geom_text(data = zerobar.df, aes(x = x, y = height, label = nlab),
                vjust = 1.4, size = list.size - 0.2, colour = "grey20", inherit.aes = FALSE) +
      expand_limits(x = zero.x - bw)
  } else if (include.zeros) {
    p <- p + expand_limits(x = zero.x - 0.1)
  }

  if (show.stats) {
    p <- p +
      geom_vline(data = stat.df, aes(xintercept = threshold),
                 linetype = "dashed", colour = "grey30", linewidth = 0.4) +
      geom_text(data = stat.df, aes(x = -Inf, y = Inf, label = label),
                hjust = -0.04, vjust = 1.3, size = stat.size, colour = "grey15",
                inherit.aes = FALSE, lineheight = 0.95)
  }

  if (label.mode %in% c("list", "both")) {
    p <- p +
      geom_text(data = header.df, aes(x = Inf, y = Inf, label = htext),
                hjust = 1.03, vjust = 1.2, size = list.size, fontface = "bold",
                colour = "grey15", inherit.aes = FALSE) +
      geom_text(data = entry.df, aes(x = Inf, y = Inf, label = entry),
                hjust = 1.03, vjust = entry.df$vpos, size = list.size,
                fontface = entry.df$fontface, colour = entry.df$col, inherit.aes = FALSE)
  }

  if (label.mode %in% c("points", "both") && requireNamespace("ggrepel", quietly = TRUE)) {
    is.bold   <- hits.show$x <= qthr.map[as.character(hits.show$Replicate)]
    is.marker <- hits.show$sgID %in% label.markers
    lab.df    <- hits.show[is.bold | is.marker, ]
    if (nrow(lab.df)) {
      lab.df$bold   <- lab.df$x <= qthr.map[as.character(lab.df$Replicate)]
      lab.df$marker <- lab.df$sgID %in% label.markers
      lab.df$face   <- ifelse(lab.df$bold, "bold", "plain")
      lab.df$col    <- ifelse(lab.df$bold & lab.df$marker,
                              grp.cols[as.character(lab.df$Group)], "black")
      p <- p + ggrepel::geom_text_repel(
        data = lab.df, aes(x = x, y = ytop, label = sgID),
        fontface = lab.df$face, colour = lab.df$col, size = list.size,
        inherit.aes = FALSE, show.legend = FALSE, max.overlaps = Inf,
        min.segment.length = 0, box.padding = 0.3, direction = "both", nudge_y = 0)
    }
  }

  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }
  return(p)
}

# ---- plot_enrichment_scatter -------------------------------------------------

#' Enrichment scatter comparing two positive replicates
#'
#' Computes a per-sgRNA enrichment score (positive cell number / unsorted cell
#' number) independently for each of two replicates, then draws a scatter plot
#' comparing the two scores. The top guides by Euclidean distance from the
#' origin are labelled by gene name.
#'
#' @inheritParams plot_neg_correlation
#' @param pos.samples Character vector of exactly two positive sample names
#'   (x-axis and y-axis).
#' @param unsorted.sample Name of the pre-sorted / unsorted denominator sample.
#' @param sg.pattern Regex stripped from sgID to extract the gene name.
#' @param n.label Number of top guides (by combined score) to label.
#' @param highlight.genes Optional character vector of gene names to label in a
#'   distinct colour, regardless of their rank. Use to cross-reference top
#'   genes from another enrichment scatter.
#' @param highlight.colour Colour for highlighted genes (default `"firebrick"`).
#' @param highlight.label Legend label for highlighted genes.
#' @param log.score If `TRUE`, apply a log transform to the enrichment scores
#'   before plotting and regression; guides with score <= 0 are dropped.
#'   If `FALSE`, plot raw ratios.
#' @param log.base Base of the log transform when `log.score = TRUE` (default
#'   10). Common choices: 10, 2, or `exp(1)` for natural log.
#' @param label.by What to use for point labels: `"gene"` (default) uses the
#'   bare gene name, `"sgID"` uses the full sgRNA identifier abbreviated from
#'   `GENE_sg5` to `GENE-5`.
#' @param min.unsorted Minimum unsorted cell number to include a guide.
#' @param label.size Size of gene-name labels.
#' @param point.size,point.alpha Aesthetic parameters.
#' @param plot.title Plot title; auto-generated if `NULL`.
#' @param save.plot If `TRUE`, save to `output.file`.
#' @param output.file File path for saved plot.
#' @param plot.width,plot.height Dimensions in inches for saved plot.
#' @return A list with `plot` (ggplot object) and `scores` (data.frame of all
#'   enrichment scores).
#' @export
plot_enrichment_scatter <- function(
    sample.info,
    pos.samples    = c("pos-1", "pos-2"),
    unsorted.sample = "unsorted-1",
    id.column      = "sgID",
    sample.column  = "Sample_ID",
    value.column   = "cell_num",
    sg.pattern     = "_sg[0-9]+$",
    n.label        = 20,
    highlight.genes = NULL,
    highlight.colour = "firebrick",
    highlight.label  = "cross-referenced",
    label.by       = c("gene", "sgID"),
    log.score      = TRUE,
    log.base       = 10,
    min.unsorted   = 1,
    label.size     = 3,
    point.size     = 0.8,
    point.alpha    = 0.35,
    palette.name   = "Set1",
    plot.title     = NULL,
    save.plot      = FALSE,
    output.file    = "enrichment_scatter.pdf",
    plot.width     = 8,
    plot.height    = 7
) {
  require(ggplot2)

  if (length(pos.samples) != 2)
    stop("pos.samples must be a character vector of exactly 2 sample names.")

  ## Read and validate
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- as.data.frame(sample.info)
  }
  needed <- c(unsorted.sample, pos.samples)
  missing.samples <- setdiff(needed, unique(df[[sample.column]]))
  if (length(missing.samples) > 0)
    stop("Sample(s) not found in '", sample.column, "': ",
         paste(missing.samples, collapse = ", "))

  ## Unsorted denominator
  uns <- df[df[[sample.column]] == unsorted.sample & df[[value.column]] >= min.unsorted,
            c(id.column, value.column)]
  names(uns) <- c("sgID", "unsorted")

  ## Enrichment per replicate
  enrich <- lapply(pos.samples, function(ps) {
    sub <- df[df[[sample.column]] == ps & df[[value.column]] > 0,
              c(id.column, value.column)]
    names(sub) <- c("sgID", "pos")
    m <- merge(sub, uns, by = "sgID")
    m$score <- m$pos / m$unsorted
    m[, c("sgID", "score")]
  })
  names(enrich) <- pos.samples

  ## Merge replicates
  merged <- merge(enrich[[1]], enrich[[2]], by = "sgID", suffixes = paste0(".", pos.samples))
  score.cols <- paste0("score.", pos.samples)
  names(merged)[names(merged) == score.cols[1]] <- "score.x"
  names(merged)[names(merged) == score.cols[2]] <- "score.y"

  ## Gene names and label column
  merged$gene <- sub(sg.pattern, "", merged$sgID)
  label.by <- match.arg(label.by)
  if (label.by == "sgID") {
    merged$label <- sub("_sg([0-9]+)$", "-\\1", merged$sgID)
  } else {
    merged$label <- merged$gene
  }

  ## Transform
  if (log.score) {
    merged <- merged[merged$score.x > 0 & merged$score.y > 0, ]
    log.fn <- function(x) log(x, base = log.base)
    merged$plot.x <- log.fn(merged$score.x)
    merged$plot.y <- log.fn(merged$score.y)
    base.lab <- if (log.base == exp(1)) "ln" else sprintf("log%g", log.base)
    x.lab <- sprintf("%s(%s / %s)", base.lab, pos.samples[1], unsorted.sample)
    y.lab <- sprintf("%s(%s / %s)", base.lab, pos.samples[2], unsorted.sample)
  } else {
    merged$plot.x <- merged$score.x
    merged$plot.y <- merged$score.y
    x.lab <- sprintf("%s / %s", pos.samples[1], unsorted.sample)
    y.lab <- sprintf("%s / %s", pos.samples[2], unsorted.sample)
  }

  ## Top N by combined score (highest enrichment across both replicates)
  merged$combined <- merged$plot.x + merged$plot.y
  merged <- merged[order(-merged$combined), ]
  top <- merged[seq_len(min(n.label, nrow(merged))), ]

  ## Report
  message(sprintf("Enrichment score: %s cell_num / %s cell_num%s",
                  "positive", unsorted.sample, if (log.score) sprintf(" (%s)", base.lab) else ""))
  message(sprintf("  %d sgRNAs with scores in both replicates", nrow(merged)))
  message(sprintf("  Top %d by distance from origin:", nrow(top)))
  print(top[, c("sgID", "gene", "score.x", "score.y", "combined")], row.names = FALSE)

  ## Title
  if (is.null(plot.title))
    plot.title <- sprintf("Enrichment: %s vs %s (normalized to %s)",
                          pos.samples[1], pos.samples[2], unsorted.sample)

  ## Highlight genes (cross-referenced from another plot)
  hl.df <- NULL
  if (!is.null(highlight.genes)) {
    hl.df <- merged[merged$gene %in% highlight.genes, ]
    hl.df <- hl.df[!duplicated(hl.df$gene), ]
  }

  ## Build plot
  top.col <- get_palette(1, palette.name)
  p <- ggplot(merged, aes(x = plot.x, y = plot.y)) +
    geom_point(size = point.size, alpha = point.alpha, colour = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_point(data = top, aes(x = plot.x, y = plot.y),
               colour = top.col, size = point.size + 1) +
    labs(x = x.lab, y = y.lab, title = plot.title) +
    theme_bw(base_size = 12)

  add_labels <- function(p, d, col, sz) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(
        data = d, aes(x = plot.x, y = plot.y, label = label),
        size = sz, colour = col,
        max.overlaps = Inf, min.segment.length = 0,
        segment.size = 0.2, box.padding = 0.35)
    } else {
      p + geom_text(data = d, aes(x = plot.x, y = plot.y, label = label),
                    size = sz, vjust = -0.6, colour = col, check_overlap = TRUE)
    }
  }
  p <- add_labels(p, top, top.col, label.size)

  if (!is.null(hl.df) && nrow(hl.df) > 0) {
    p <- p +
      geom_point(data = hl.df, aes(x = plot.x, y = plot.y),
                 colour = highlight.colour, size = point.size + 1.2, shape = 17)
    p <- add_labels(p, hl.df[!(hl.df$gene %in% top$gene), ], highlight.colour, label.size)
  }

  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  invisible(list(plot = p, scores = merged))
}

# ---- plot_gene_set_bars ------------------------------------------------------

#' Horizontal bar chart of top gene sets
#'
#' Visualises gene-set membership results from an MSigDB lookup.
#' Input is a data.frame with columns `gene_set`, `n_genes`, and `genes`
#' (as returned by a `gene_set_lookup()` helper).  Bars are ordered by
#' `n_genes` and annotated with the matching gene names.
#'
#' @param set_counts A data.frame with columns \code{gene_set} (character),
#'   \code{n_genes} (integer), and \code{genes} (comma-separated character
#'   string of matching gene symbols).
#' @param n_top Maximum number of gene sets to show (default 15).
#' @param min_genes Minimum number of matching genes for a set to be included
#'   (default 2).
#' @param title Plot title.
#' @param bar.colour Fill colour for the bars.
#' @param label.size Font size for gene-name annotations on bars.
#' @param save.plot If \code{TRUE}, save to \code{output.file}.
#' @param output.file File path for saved plot.
#' @param plot.width,plot.height Dimensions in inches for saved plot.
#' @return A ggplot object (invisibly).
#' @export
plot_gene_set_bars <- function(
    set_counts,
    n_top         = 15,
    min_genes     = 2,
    title         = "Top gene sets",
    bar.colour    = "steelblue",
    label.size    = 3,
    save.plot     = FALSE,
    output.file   = "gene_set_bars.pdf",
    plot.width    = 10,
    plot.height   = 6
) {
  if (is.null(set_counts) || nrow(set_counts) == 0) {
    message("No gene set data to plot.")
    return(invisible(NULL))
  }

  ## Filter and trim
  df <- set_counts[set_counts$n_genes >= min_genes, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No gene sets meet the min_genes threshold (", min_genes, ").")
    return(invisible(NULL))
  }
  df <- head(df[order(-df$n_genes), ], n_top)

  ## Wrap long set names for readability
  df$gene_set_wrap <- sapply(df$gene_set, function(s) {
    paste(strwrap(tolower(s), width = 40), collapse = "\n")
  })

  ## Order factor by n_genes (lowest at top so highest prints at top of horiz bar)
  df$gene_set_wrap <- factor(df$gene_set_wrap,
                             levels = rev(df$gene_set_wrap))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_genes, y = gene_set_wrap)) +
    ggplot2::geom_col(fill = bar.colour, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = genes),
                       hjust = -0.05, size = label.size, colour = "grey20") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.55)),
      breaks = function(lim) seq(0, floor(lim[2]), by = 1)
    ) +
    ggplot2::labs(x = "Matching genes", y = NULL, title = title) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  if (save.plot) {
    ggplot2::ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  invisible(p)
}

# ---- plot_enrichment_barcode -------------------------------------------------

#' Per-gene barcode plot of sgRNA enrichment scores
#'
#' Shows the top genes ranked by median enrichment score, with each sgRNA
#' drawn as an individual bar so you can assess whether enrichment is
#' consistent across multiple guides or driven by a single outlier.
#' Enrichment is computed as sample cell number / unsorted cell number,
#' averaged across the supplied replicates.
#'
#' @inheritParams plot_enrichment_scatter
#' @param n_genes Number of top genes to display (default 20).
#' @param bar.width Width of individual sgRNA bars (default 0.7).
#' @param bar.colour Fill colour for the bars.
#' @param median.colour Colour of the median indicator line per gene.
#' @param title Plot title; auto-generated if `NULL`.
#' @param save.plot If `TRUE`, save to `output.file`.
#' @param output.file File path for saved plot.
#' @param plot.width,plot.height Dimensions in inches for saved plot.
#' @return A list with `plot` (ggplot object) and `gene_scores` (data.frame
#'   of per-gene median scores, sorted descending).
#' @export
plot_enrichment_barcode <- function(
    sample.info,
    pos.samples     = c("pos-1", "pos-2"),
    unsorted.sample = "unsorted-1",
    id.column       = "sgID",
    sample.column   = "Sample_ID",
    value.column    = "cell_num",
    sg.pattern      = "_sg[0-9]+$",
    n_genes         = 20,
    log.score       = TRUE,
    log.base        = 10,
    min.unsorted    = 1,
    bar.width       = 0.7,
    bar.colour      = "steelblue",
    median.colour   = "firebrick",
    title           = NULL,
    save.plot       = FALSE,
    output.file     = "enrichment_barcode.pdf",
    plot.width      = 10,
    plot.height     = 7
) {
  require(ggplot2)

  ## Read data
  if (is.character(sample.info)) {
    df <- read.csv(sample.info, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- as.data.frame(sample.info)
  }

  ## Unsorted denominator
  uns <- df[df[[sample.column]] == unsorted.sample & df[[value.column]] >= min.unsorted,
            c(id.column, value.column)]
  names(uns) <- c("sgID", "unsorted")

  ## Enrichment per replicate, then average across replicates per sgRNA
  enrich_list <- lapply(pos.samples, function(ps) {
    sub <- df[df[[sample.column]] == ps & df[[value.column]] > 0,
              c(id.column, value.column)]
    names(sub) <- c("sgID", "pos")
    m <- merge(sub, uns, by = "sgID")
    m$score <- m$pos / m$unsorted
    m[, c("sgID", "score")]
  })
  # Merge all replicates
  scores <- enrich_list[[1]]
  names(scores)[2] <- "score.1"
  for (i in seq_along(enrich_list)[-1]) {
    el <- enrich_list[[i]]
    names(el)[2] <- paste0("score.", i)
    scores <- merge(scores, el, by = "sgID")
  }
  score.cols <- paste0("score.", seq_along(pos.samples))
  scores$mean_score <- rowMeans(scores[, score.cols, drop = FALSE], na.rm = TRUE)

  ## Gene and sgRNA label
  scores$gene  <- sub(sg.pattern, "", scores$sgID)
  scores$sg_id <- sub("_sg([0-9]+)$", "-\\1", scores$sgID)

  ## Log transform
  if (log.score) {
    scores <- scores[scores$mean_score > 0, ]
    scores$plot_score <- log(scores$mean_score, base = log.base)
    base.lab <- if (log.base == exp(1)) "ln" else sprintf("log%g", log.base)
    y.lab <- sprintf("%s(enrichment score)", base.lab)
  } else {
    scores$plot_score <- scores$mean_score
    y.lab <- "Enrichment score"
  }

  ## Per-gene median
  gene_med <- aggregate(plot_score ~ gene, data = scores, FUN = median)
  names(gene_med)[2] <- "median_score"
  gene_med <- gene_med[order(-gene_med$median_score), ]

  ## Top N genes
  top_genes <- head(gene_med$gene, n_genes)
  plot_df <- scores[scores$gene %in% top_genes, ]
  plot_df$gene <- factor(plot_df$gene, levels = rev(top_genes))
  gene_med_top <- gene_med[gene_med$gene %in% top_genes, ]
  gene_med_top$gene <- factor(gene_med_top$gene, levels = rev(top_genes))

  ## Title
  if (is.null(title))
    title <- sprintf("Top %d genes by median enrichment (%s)",
                     length(top_genes), paste(pos.samples, collapse = " + "))

  ## Plot — horizontal bars per sgRNA, genes on y-axis
  p <- ggplot(plot_df, aes(x = plot_score, y = gene)) +
    geom_point(aes(colour = sg_id), size = 3, shape = 124, stroke = 3) +
    geom_point(data = gene_med_top, aes(x = median_score, y = gene),
               colour = median.colour, shape = 18, size = 3) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    labs(x = y.lab, y = NULL, title = title, colour = "sgRNA") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none"
    )

  ## Add sgRNA labels next to each tick mark
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = sg_id), size = 2.5, colour = "grey30",
      direction = "x", nudge_x = 0.05,
      max.overlaps = Inf, segment.size = 0.15, box.padding = 0.15
    )
  }

  if (save.plot) {
    ggsave(output.file, plot = p, width = plot.width, height = plot.height)
    message("Saved plot to: ", output.file)
  }

  invisible(list(plot = p, gene_scores = gene_med))
}
