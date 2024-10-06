library(ape)

tr <- read.nexus('socher_albanoRomancePW.mcc.tree')
write.tree(tr,'socher_albanoRomancePW.mcc.nwk')
tr <- read.nexus('socher_albanoRomanceSW.mcc.tree')
write.tree(tr,'socher_albanoRomanceSW.mcc.nwk')
tr <- read.nexus('socher_albanoRomance.mcc.tree')
write.tree(tr,'socher_albanoRomance.mcc.nwk')


