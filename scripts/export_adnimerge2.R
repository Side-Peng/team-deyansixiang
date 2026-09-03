# export_adnimerge2.R — 将 ADNIMERGE2 的 .rda 表一次性导出为 CSV（管线所需 10 张）
# 运行：Rscript scripts/export_adnimerge2.R （在 goai-open-exploration/ 下）
# 输出：data/raw/merged/<TABLE>.csv + 每张表的列名打印（用于核对 config.py 字段映射）

tables <- c(
  "PTDEMOG", "DXSUM", "CDR", "MMSE", "ADAS", "MEDHIST", "INITHEALTH",
  "UPENN_PLASMA_FUJIREBIO_QUANTERIX", "UCBERKELEY_AMY_6MM", "UCBERKELEY_TAU_6MM"
)

src <- "ADNI_MERGE/ADNIMERGE2/data"
out <- "data/raw/merged"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

for (t in tables) {
  e <- new.env()
  load(file.path(src, paste0(t, ".rda")), envir = e)
  obj <- ls(e)[1]
  d <- get(obj, envir = e)
  cat("==", t, "==", nrow(d), "rows,", ncol(d), "cols\n")
  cat("  names:", paste(names(d), collapse = ", "), "\n")
  write.csv(d, file.path(out, paste0(t, ".csv")), row.names = FALSE, na = "")
}

cat("\nDONE ->", normalizePath(out), "\n")
