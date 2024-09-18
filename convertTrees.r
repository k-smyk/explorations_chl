library(ape)

t1 <- read.nexus('albanoRomancePW.run1.t')
write.tree(t1,'albanoRomancePW.run1.tre')

t2 <- read.nexus('albanoRomancePW.run2.t')
write.tree(t2,'albanoRomancePW.run2.tre')

t1 <- read.nexus('albanoRomanceSW.run1.t')
write.tree(t1,'albanoRomanceSW.run1.tre')

t2 <- read.nexus('albanoRomanceSW.run2.t')
write.tree(t2,'albanoRomanceSW.run2.tre')
