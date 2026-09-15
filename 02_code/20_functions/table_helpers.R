# table_helpers.R
# Shared rendering for the manuscript tables. Every table is built as a data.table and then
# passed through BOTH of these, so the HTML you review and the Word file you submit come from
# one object and can't drift apart.
#
# House style: Arial, no captions (those live in the manuscript), section headers bold and
# carrying their own total, children indented, landscape because these are wide.

suppressMessages({library(flextable); library(officer); library(knitr); library(kableExtra)})

# a credible/confidence interval that excludes the null, for bolding.
# `null_at` is 0 for differences and 1 for ratios.
sig_cells <- function(dt, pattern = "CrI|CI", null_at = 0) {
  m <- matrix(FALSE, nrow(dt), ncol(dt))
  rows <- if (is.null(pattern)) seq_len(nrow(dt)) else grep(pattern, dt[[1]])
  for (j in 2:ncol(dt)) {
    if (!is.character(dt[[j]])) next          # factors and numerics carry no interval
    for (i in rows) {
    v <- dt[[j]][i]
    if (is.na(v) || !nzchar(v)) next
    g <- regmatches(v, regexec("^(-?[0-9.]+) \\((-?[0-9.]+), (-?[0-9.]+)\\)$", v))[[1]]
    if (length(g) == 4) {
      lo <- as.numeric(g[3]); hi <- as.numeric(g[4])
      m[i, j] <- (lo > null_at) || (hi < null_at)
    }
  }}
  m
}

# `span` is a named vector giving a header row ABOVE the column names, e.g.
#   c(" " = 1, "Cumulative rate ratios" = 4, "Lag-specific rate ratios" = 4)
show_tbl <- function(dt, indent_rows = integer(0), bold_rows = integer(0),
                     sig = FALSE, null_at = 0, pattern = "CrI|CI", font_size = 13,
                     span = NULL, indent2 = integer(0), col_labels = NULL, vrule = NULL,
                     title = NULL, note = NULL) {
  d <- data.table::copy(dt)
  if (sig) {
    m <- sig_cells(d, pattern, null_at)
    for (j in 2:ncol(d)) d[[j]][m[, j]] <- paste0("**", d[[j]][m[, j]], "**")
  }
  # display labels can repeat (e.g. the same lag windows under two spanning headers) even
  # though the underlying column names must stay unique
  k <- kbl(d, align = c("l", rep("r", ncol(d) - 1)), booktabs = TRUE,
           caption = title,
           col.names = if (is.null(col_labels)) names(d) else col_labels) |>
    kable_styling(bootstrap_options = "condensed", full_width = FALSE,
                  font_size = font_size, html_font = "Arial, Helvetica, sans-serif")
  if (!is.null(vrule))     k <- column_spec(k, vrule, border_right = TRUE)
  if (!is.null(span))      k <- add_header_above(k, span, bold = TRUE)
  if (length(indent_rows)) k <- add_indent(k, indent_rows)
  if (length(indent2))     k <- add_indent(k, indent2, level_of_indent = 2)
  if (length(bold_rows))   k <- row_spec(k, bold_rows, bold = TRUE)
  k <- row_spec(k, 0, bold = TRUE)
  if (!is.null(note)) k <- kableExtra::footnote(k, general = note, general_title = "", footnote_as_chunk = FALSE)
  k
}

to_word <- function(dt, name, folder, indent_rows = integer(0), bold_rows = integer(0),
                    rule_after = integer(0), sig = FALSE, null_at = 0,
                    pattern = "CrI|CI", font_size = 10, max_width = 9.4,
                    span = NULL, indent2 = integer(0), col_labels = NULL, vrule = NULL,
                    title = NULL, note = NULL) {
  ft <- flextable(as.data.frame(dt))
  if (!is.null(col_labels))
    ft <- set_header_labels(ft, values = setNames(as.list(col_labels), names(dt)))
  ft <- ft |> fontsize(size = font_size, part = "all") |> font(fontname = "Arial", part = "all") |>
    align(j = 1, align = "left", part = "all") |>
    align(j = 2:ncol(dt), align = "right", part = "all") |>
    bold(part = "header") |>
    padding(padding.top = 2, padding.bottom = 2, padding.right = 4, part = "all") |>
    theme_booktabs() |>
    hline_top(border = fp_border(color = "black", width = 1.2), part = "header") |>
    hline_bottom(border = fp_border(color = "black", width = 1.2), part = "body")
  if (!is.null(span)) {
    ft <- add_header_row(ft, values = names(span), colwidths = as.integer(span), top = TRUE) |>
      bold(part = "header") |> align(align = "center", part = "header") |>
      hline_top(border = fp_border(color = "black", width = 1.2), part = "header")
    for (k in which(names(span) != " "))
      ft <- hline(ft, i = 1, j = (cumsum(c(0, span))[k] + 1):cumsum(span)[k],
                  border = fp_border(color = "black", width = 0.8), part = "header")
  }
  if (length(bold_rows))   ft <- bold(ft, i = bold_rows, j = 1)
  if (length(indent_rows)) ft <- padding(ft, i = indent_rows, j = 1, padding.left = 12)
  if (length(indent2))     ft <- padding(ft, i = indent2, j = 1, padding.left = 26)
  if (length(rule_after))  ft <- hline(ft, i = rule_after,
                                       border = fp_border(color = "grey70", width = 0.5))
  if (!is.null(vrule))     ft <- vline(ft, j = vrule,
                                       border = fp_border(color = "grey40", width = 1),
                                       part = "all")
  if (sig) {
    m <- sig_cells(dt, pattern, null_at)
    for (j in 2:ncol(dt)) { i <- which(m[, j]); if (length(i)) ft <- bold(ft, i = i, j = j) }
  }
  # footnote goes under the table; the title is a paragraph above it, per journal style
  if (!is.null(note))
    ft <- add_footer_lines(ft, note) |>
      fontsize(size = max(7, font_size - 1), part = "footer") |>
      font(fontname = "Arial", part = "footer") |>
      align(align = "left", part = "footer")
  ft <- autofit(ft) |> fit_to_width(max_width = max_width)
  doc <- read_docx()
  if (!is.null(title))
    doc <- body_add_fpar(doc, fpar(ftext(title, fp_text(font.family = "Arial", font.size = 10,
                                                        bold = TRUE))))
  # stash the table AND how it should be displayed, so one brief can re-render every exhibit
  # exactly as it was exported rather than cross-referencing other documents
  try(saveRDS(list(dt = dt, indent_rows = indent_rows, indent2 = indent2, bold_rows = bold_rows,
                   span = span, col_labels = col_labels, title = title, note = note,
                   sig = sig, null_at = null_at, pattern = pattern),
              file.path(dirname(folder), "table_specs", paste0(name, ".RDS"))), silent = TRUE)
  doc |> body_add_flextable(ft) |> body_end_section_landscape() |>
    print(target = file.path(folder, paste0(name, ".docx")))
  invisible(dt)
}

# Build a sectioned table: a bold heading row per group, an indented sub-heading per item, and
# indented data rows beneath. `get` returns a character vector of length(cols) or NULL to skip.
# Returns the table plus the row indices each level occupies, so the caller can pass them
# straight to show_tbl()/to_word() rather than recomputing with which().
build_sections <- function(groups, items, leaves, get, cols, item_vals = NULL) {
  rows <- list(); lv <- integer(0)
  add <- function(lab, vals, level) {
    rows[[length(rows) + 1L]] <<- as.list(c(setNames(lab, " "),
      setNames(if (is.null(vals)) rep("", length(cols)) else as.character(vals), cols)))
    lv[length(lv) + 1L] <<- level
  }
  for (g in groups) {
    add(g, NULL, 1L)
    for (it in items) {
      leaf_vals <- lapply(leaves, function(lf) get(g, it, lf))
      if (all(vapply(leaf_vals, is.null, TRUE))) next
      add(it, if (is.null(item_vals)) NULL else item_vals(g, it), 2L)
      for (k in seq_along(leaves)) if (!is.null(leaf_vals[[k]])) add(leaves[k], leaf_vals[[k]], 3L)
    }
  }
  list(dt = data.table::rbindlist(rows, fill = TRUE),
       g1 = which(lv == 1L), g2 = which(lv == 2L), g3 = which(lv == 3L))
}

# Re-render a table exactly as to_word() exported it. to_word() stashes the data plus every
# display argument in 04_tables/table_specs/, so the manuscript brief can show the real exhibit
# instead of pointing at another html file. Returns NULL if the spec hasn't been built yet.
show_spec <- function(name, specs_dir = NULL, font_size = 12) {
  p <- file.path(specs_dir %||% file.path(tables_folder, "table_specs"), paste0(name, ".RDS"))
  if (!file.exists(p)) return(invisible(NULL))
  s <- readRDS(p)
  show_tbl(s$dt, indent_rows = s$indent_rows %||% integer(0),
           indent2 = s$indent2 %||% integer(0), bold_rows = s$bold_rows %||% integer(0),
           span = s$span, col_labels = s$col_labels, title = s$title, note = s$note,
           sig = isTRUE(s$sig), null_at = s$null_at %||% 0,
           pattern = s$pattern %||% "CrI|CI", font_size = font_size)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
