library(ggplot2)
library(dplyr)
library(ggrepel)
args <- commandArgs(trailingOnly = TRUE)
UPpeak <- read.csv(args[1], header = TRUE)
Downpeak <- read.csv(args[2], header = TRUE)

P <- data.frame(Group = character(), Value = numeric(), stringsAsFactors = FALSE)

keywords <- c("Promoter", "5' UTR", "Distal Intergenic", "Intron", "Exon", "3' UTR")

for(keyword in keywords) {
   count <- sum(grepl(keyword, UPpeak$annotation, ignore.case = TRUE))
   group_name <- paste("Increased", gsub("'", "", keyword), sep = "_")
   P <- rbind(P, data.frame(Group = group_name, Value = count, stringsAsFactors = FALSE))
}

for(keyword in keywords) {
   count <- sum(grepl(keyword, Downpeak$annotation, ignore.case = TRUE))
   group_name <- paste("Decreased", gsub("'", "", keyword), sep = "_")
   P <- rbind(P, data.frame(Group = group_name, Value = -count, stringsAsFactors = FALSE))
}

P$Value <- as.numeric(P$Value)
pdf(args[3],width=10,height=7)
plot <- ggplot(P, aes(x = reorder(Group, Value), y = Value, fill = Value)) +
   geom_bar(stat = "identity", show.legend = FALSE) +
   geom_text(aes(label = round(Value, 1),
                 hjust = ifelse(Value < 0, 1.2, -0.2),
                 vjust = 0.5),
             size = 3) +
   xlab("Types_of_peaks") +
   ylab("Number_of_peaks") +
   coord_flip() + 
   scale_y_continuous(limits = c(min(P$Value) * 1.1, max(P$Value) * 1.1))+
   scale_fill_gradient2(low = "#6495ED",
                        mid = "aliceblue",
                        high = "#FF0000")+
   theme_minimal() +  
   theme(
      panel.grid = element_blank(),  
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.ticks = element_line(size = 0.5), 
      axis.text = element_text(size = 12),  
      axis.title = element_text(size = 14) 
   ) 
print(plot)
dev.off()

keywords1 <- c("Promoter", "5 UTR", "Distal Intergenic", "Intron", "Exon", "3 UTR")
data_chars <- data.frame(Group = keywords1, 
                 AllPeaks = sapply(keywords1, function(keyword) {
                    sum(abs(P$Value[grepl(keyword, P$Group, ignore.case = TRUE)]))
                 }),
                 stringsAsFactors = FALSE)
pdf(args[4],width=5,height=4)
data_chars %>% ggplot(aes(x = reorder(Group,-AllPeaks), y = AllPeaks)) + geom_bar(stat = "identity",width = 0.5,fill = "#0099FF")+theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),axis.line = element_line(color = "black",linewidth = 1),panel.background = element_blank(),                        panel.grid = element_blank())

dev.off()
