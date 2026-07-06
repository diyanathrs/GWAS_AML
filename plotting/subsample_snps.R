subsample_gwas <- function(df,
                           pcol = "p",
                           chrcol = "chr",
                           poscol = "pos",
                           keep_threshold = 1e-4,
                           max_per_chr = 5000,
                           random_seed = 1234) {
  set.seed(random_seed)
  
  # Ensure required columns exist
  if (!all(c(pcol, chrcol, poscol) %in% names(df))) {
    stop("Missing required columns in df.")
  }
  
  # Compute -log10(p)
  df$logp <- -log10(df[[pcol]])
  
  # 1. Keep all significant SNPs
  df_keep <- df[df[[pcol]] <= keep_threshold, ]
  
  # 2. Subsample the rest *within each chromosome*
  df_sub <- df[df[[pcol]] > keep_threshold, ]
  
  df_subsampled <- do.call(rbind, lapply(split(df_sub, df_sub[[chrcol]]), function(x) {
    n <- nrow(x)
    if (n <= max_per_chr) {
      return(x)
    } else {
      idx <- sample.int(n, max_per_chr)
      return(x[idx, ])
    }
  }))
  
  # 3. Combine and sort by chromosome & position
  result <- rbind(df_keep, df_subsampled)
  result <- result[order(result[[chrcol]], result[[poscol]]), ]
  
  rownames(result) <- NULL
  return(result)
}