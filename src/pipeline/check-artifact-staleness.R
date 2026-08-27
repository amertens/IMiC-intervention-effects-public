# =============================================================================
# check-artifact-staleness.R
#
# Flag any manuscript artifact that is OLDER than something it is derived from.
# This is the check that would have caught the Fig. S10 failure -- a figure dated
# 2026-08-12 plotting a CSV regenerated on 2026-08-14 -- before it reached the
# supplement, and it catches the same class of drift for every other exhibit.
#
# Reads results/ARTIFACT_MANIFEST.csv (run src/pipeline/build-artifact-manifest.R
# first; this script will build it if absent). Three checks per artifact:
#   STALE-DATA  : the artifact is older than one of its declared `depends` inputs
#   STALE-CODE  : the artifact is older than the generator script that writes it
#   NO-GENERATOR: nothing in the repo produces it (a static/imported binary)
#   MISSING     : the declared path does not exist
#
# Exit status is 1 if any STALE-* or MISSING row is found, so this can gate a
# release build; NO-GENERATOR rows are reported but do not fail the run (Figs S1,
# S2, S5, S6 are deliberately static).
#
# Run from repo root: Rscript src/pipeline/check-artifact-staleness.R
# =============================================================================
suppressMessages({ library(data.table) })
root <- paste0(here::here(), "/")
MAN  <- paste0(root, "results/ARTIFACT_MANIFEST.csv")
OUT  <- paste0(root, "results/ARTIFACT_STALENESS_REPORT.csv")

# A generated file can legitimately be a few seconds older than an input touched in
# the same pipeline run (they are written sequentially), so only differences beyond
# this margin count as stale.
TOLERANCE_SECS <- 120

if (!file.exists(MAN)) {
  message("manifest absent; building it first")
  source(paste0(root, "src/pipeline/build-artifact-manifest.R"))
}
man <- fread(MAN)

mtime_of <- function(p) {
  # a path field may list several files joined by ";" (e.g. Table S3 = combined +
  # stratified CSV). Treat the artifact as present only if ALL parts exist, and
  # take the OLDEST part's mtime (conservative: any stale component => stale).
  parts <- trimws(unlist(strsplit(p, ";", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  f <- file.path(root, parts)
  if (!length(f) || !all(file.exists(f))) return(NA_real_)
  min(as.numeric(file.info(f)$mtime))
}

rows <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  a <- man[i]
  at <- mtime_of(a$path)
  if (is.na(at))
    return(data.table(status = "MISSING", exhibit = a$exhibit, artifact = a$path,
                      input = "", artifact_mtime = "", input_mtime = "", lag_hours = NA_real_))
  out <- list()
  # code staleness
  if (nzchar(a$generator)) {
    gt <- mtime_of(a$generator)
    if (!is.na(gt) && gt - at > TOLERANCE_SECS)
      out[[length(out) + 1]] <- data.table(
        status = "STALE-CODE", exhibit = a$exhibit, artifact = a$path, input = a$generator,
        artifact_mtime = format(as.POSIXct(at, origin = "1970-01-01"), "%Y-%m-%d %H:%M"),
        input_mtime    = format(as.POSIXct(gt, origin = "1970-01-01"), "%Y-%m-%d %H:%M"),
        lag_hours = round((gt - at) / 3600, 1))
  } else {
    out[[length(out) + 1]] <- data.table(
      status = "NO-GENERATOR", exhibit = a$exhibit, artifact = a$path, input = "",
      artifact_mtime = format(as.POSIXct(at, origin = "1970-01-01"), "%Y-%m-%d %H:%M"),
      input_mtime = "", lag_hours = NA_real_)
  }
  # data staleness
  deps <- trimws(unlist(strsplit(a$depends, ";", fixed = TRUE)))
  for (d in deps[nzchar(deps)]) {
    dt <- mtime_of(d)
    if (is.na(dt)) {
      out[[length(out) + 1]] <- data.table(
        status = "MISSING", exhibit = a$exhibit, artifact = a$path, input = d,
        artifact_mtime = "", input_mtime = "", lag_hours = NA_real_)
    } else if (dt - at > TOLERANCE_SECS) {
      out[[length(out) + 1]] <- data.table(
        status = "STALE-DATA", exhibit = a$exhibit, artifact = a$path, input = d,
        artifact_mtime = format(as.POSIXct(at, origin = "1970-01-01"), "%Y-%m-%d %H:%M"),
        input_mtime    = format(as.POSIXct(dt, origin = "1970-01-01"), "%Y-%m-%d %H:%M"),
        lag_hours = round((dt - at) / 3600, 1))
    }
  }
  if (!length(out)) return(NULL)
  rbindlist(out)
}))

if (!nrow(rows)) rows <- data.table(status = character(), exhibit = character(),
  artifact = character(), input = character(), artifact_mtime = character(),
  input_mtime = character(), lag_hours = numeric())
setorder(rows, status, exhibit)
fwrite(rows, OUT)

fail <- rows[status %in% c("STALE-DATA", "STALE-CODE", "MISSING")]
cat("artifact staleness report ->", OUT, "\n")
cat("  checked:", nrow(man), "artifacts\n")
for (s in c("STALE-DATA", "STALE-CODE", "MISSING", "NO-GENERATOR"))
  cat(sprintf("  %-13s %d\n", s, sum(rows$status == s)))
if (nrow(fail)) {
  cat("\n--- needs a rebuild ---\n")
  print(fail[, .(status, exhibit, artifact, input, artifact_mtime, input_mtime, lag_hours)],
        row.names = FALSE)
}
if (nrow(rows[status == "NO-GENERATOR"])) {
  cat("\n--- no generator in repo (static/imported; not a failure) ---\n")
  print(rows[status == "NO-GENERATOR", .(exhibit, artifact)], row.names = FALSE)
}
quit(status = if (nrow(fail)) 1L else 0L)
