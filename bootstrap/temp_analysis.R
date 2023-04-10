library(ggplot2)

x <- read.csv("bootstrap/input/OX_WT.tsv", sep="\t")
x <- x[(x$Read=="R1") & (x$Name=="without_AI") & (x$Breaks=="2dsb") & (x$Construct %in% c("Sense", "BranchD")), ]

y_a_a <- x[(x$Cell_line=="WT") & (x$Construct=="Sense"), ]
y_a_b <- x[(x$Cell_line=="WT") & (x$Construct=="BranchD"), ]
y_b_a <- x[(x$Cell_line=="OX") & (x$Construct=="Sense"), ]
y_b_b <- x[(x$Cell_line=="OX") & (x$Construct=="BranchD"), ]

data = data.frame(
  x=factor(c("WT", "OX"), c("WT", "OX")),
  y=c(mean(y_a_a$Frequency/y_a_b$Frequency), mean(y_b_a$Frequency/y_b_b$Frequency)),
  sd=c(sd(y_a_a$Frequency/y_a_b$Frequency), sd(y_b_a$Frequency/y_b_b$Frequency))
) |>
ggplot(data = _, mapping = aes(x, y, ymin = y-sd, ymax = y + sd)) +
geom_bar(stat = "identity") +
geom_errorbar()
ggsave("ratio_then_mean.png", height=3, width=3)

data = data.frame(
  x=factor(c("WT", "OX"), c("WT", "OX")),
  y=c(mean(y_a_a$Frequency)/mean(y_a_b$Frequency), mean(y_b_a$Frequency)/mean(y_b_b$Frequency))
) |>
ggplot(data = _, mapping = aes(x, y)) +
geom_bar(stat = "identity")
ggsave("mean_then_ratio.png", height=3, width=3)