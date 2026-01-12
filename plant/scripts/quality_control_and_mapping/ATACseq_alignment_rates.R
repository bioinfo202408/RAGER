library(tidyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

file_paths <- strsplit(args[1], " ")[[1]]

output_pdf <- args[2]

samples <- sapply(file_paths, function(x) {
   parts <- strsplit(x, "/")[[1]]
   basename(parts[length(parts) - 1])
})

p <- data.frame(Sample = character(), UnAR = numeric(), UniqueAR = numeric(), MultiAR = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(file_paths)) {
   file_path <- file_paths[i]
   
   lines <- readLines(file_path)
   
   alignment_lines <- grep("aligned concordantly", lines, value = TRUE)
   UnAR <- as.numeric(sub(".*\\((.*?)%\\).*", "\\1", alignment_lines[1]))
   UniqueAR <- as.numeric(sub(".*\\((.*?)%\\).*", "\\1", alignment_lines[2]))
   MultiAR <- as.numeric(sub(".*\\((.*?)%\\).*", "\\1", alignment_lines[3]))
   
   temp_df <- data.frame(Sample = samples[i], UnAR = UnAR, UniqueAR = UniqueAR, MultiAR = MultiAR)
   p <- rbind(p, temp_df)
}

p_long <- pivot_longer(p,
                       cols = c("UnAR", "UniqueAR", "MultiAR"),
                       names_to = "Category",
                       values_to = "Rate")

pdf(output_pdf, width = 6, height = 6)
ggplot(p_long, aes(x = Sample, y = Rate, fill = Category)) +
   geom_bar(stat = "identity", position = "stack") +
   labs(x = "Sample", y = "Rate") +
   scale_fill_manual(values = c("steelblue", "orange", "#A5CCE5")) +
   theme_minimal() +
   ggtitle("ATACseq Alignment Rates")
dev.off()
