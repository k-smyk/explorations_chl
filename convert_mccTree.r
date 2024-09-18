library(ape)

tr <- read.nexus('albanoRomancePW.mcc.tree')
write.tree(tr,'albanoRomancePW.mcc.nwk')
tr <- read.nexus('albanoRomanceSW.mcc.tree')
write.tree(tr,'albanoRomanceSW.mcc.nwk')


