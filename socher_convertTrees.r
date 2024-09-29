library(ape)

t1 <- read.nexus('socher_albanoRomancePW.run1.t')
write.tree(t1,'socher_albanoRomancePW.run1.tre')

t2 <- read.nexus('socher_albanoRomancePW.run2.t')
write.tree(t2,'socher_albanoRomancePW.run2.tre')

t1 <- read.nexus('socher_albanoRomanceSW.run1.t')
write.tree(t1,'socher_albanoRomanceSW.run1.tre')

t2 <- read.nexus('socher_albanoRomanceSW.run2.t')
write.tree(t2,'socher_albanoRomanceSW.run2.tre')
