x <- read.csv("C:\\Users\\tchan\\Code\\DSB_Project\\RNA-mediated_DSB_repair\\test_ratios\\input\\mmej.tsv", sep="\t")

x[(x$Read == "R1") & (x$Name == "EE1") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R1") & (x$Name == "EE2") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R1") & (x$Name == "EE3") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R1") & (x$Name == "EE9") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R1") & (x$Name == "EI2") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R2") & (x$Name == "EE1R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R2") & (x$Name == "EI28R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD")), ]
x[(x$Read == "R2") & (x$Name == "EI30R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD")), ]


x[
  !(
    ((x$Read == "R1") & (x$Name == "EE1") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R1") & (x$Name == "EE2") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R1") & (x$Name == "EE3") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R1") & (x$Name == "EE9") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R1") & (x$Name == "EI2") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R2") & (x$Name == "EE1R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R2") & (x$Name == "EI28R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R2") & (x$Name == "EI30R") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD")))
  ) &
  (
    ((x$Read == "R1") & (x$Breaks == "sgA") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R2") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD"))) |
    ((x$Read == "R2") & (x$Breaks == "2dsb") & (x$Construct %in% c("Sense", "BranchD")))
  ) &
  !startsWith(x$Name, "EB"),
] -> y

y[(y[, "Frequency"] == 0), ]

# MMMEJ w/ at least 1 out of 16 zero frequency
# sgA, R1, EE1,
# sgA, R1, EE2,
# sgA, R1, EE3,
# sgA, R1, EE9,
# sgA, R1, EI2,
# sgB, R2, EE1R,
# sgB, R2, EI28R,
# sgB, R2, EI30R,

x <- read.csv("C:\\Users\\tchan\\Code\\DSB_Project\\RNA-mediated_DSB_repair\\test_ratios\\input\\R-TDR.tsv", sep="\t")
x[(x$Read == "R2") & (x$Name == "R-TDR") & (x$Breaks == "sgB") & (x$Construct %in% c("Sense", "BranchD")), ]
