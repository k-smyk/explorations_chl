import pandas as pd
from numpy import *


def nexCharOutput(chMtx,names,outfile,datatype='STANDARD'):
    f = open(outfile,'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax='+str(len(chMtx))+' NCHAR='+str(len(chMtx.T))+';\n')
    f.write('FORMAT DATATYPE='+datatype+' GAP=? MISSING=- interleave=yes;\n')
    f.write('MATRIX\n\n')
    txLgth = max(map(len,names))
    for i in range(len(chMtx)):
        f.write(names[i].ljust(txLgth+2))
        for ch in chMtx[i]:
            if ch==-1: ch='-'
            else:
                ch = str(ch)
            f.write(ch)
        f.write('\n')
    f.write('\n;\n\nEND;\n')
    f.close()


scPW = pd.read_table('pairwise_alignment/albanoRomanceSC-pw.nex',
                   skiprows=7,
                   skipfooter=4,engine='python',
                   index_col=0,
                   sep='\s+',header=None)

scSW = pd.read_table('pairwise_alignment/albanoRomanceSC-sw.nex',
                   skiprows=7,
                   skipfooter=4,engine='python',
                   index_col=0,
                   sep='\s+',header=None)


ccPW = pd.read_table('albanoRomanceCC_PW.nex',
                   skiprows=7,
                   skipfooter=4,engine='python',
                   index_col=0,
                   sep='\s+',header=None)

ccSW = pd.read_table('albanoRomanceCC_SW.nex',
                   skiprows=7,
                   skipfooter=4,engine='python',
                   index_col=0,
                   sep='\s+',header=None)

ccPW = ccPW.loc[scPW.index]
ccSW = ccSW.loc[scSW.index]

scCC_PW = pd.DataFrame([list(scPW[1][l]+ccPW[1][l]) for l in scPW.index],
                    index = scPW.index)
scCC_SW = pd.DataFrame([list(scSW[1][l]+ccSW[1][l]) for l in scSW.index],
                    index = scSW.index)

nexCharOutput(scCC_PW.values,scCC_PW.index,'albanoRomance_sc_cc_PW.nex',datatype='restriction')
nexCharOutput(scCC_SW.values,scCC_SW.index,'albanoRomance_sc_cc_SW.nex',datatype='restriction')

nPW = len(scPW.values[0][0])
nSW = len(scSW.values[0][0])
mPW = len(ccPW.values[0][0])
mSW = len(ccSW.values[0][0])

mbCommandsPW = """#NEXUS
begin MrBayes;
      execute albanoRomance_sc_cc_PW.nex;
      charset sc = 1-"""+str(nPW)+""";
      charset cc = """+str(nPW+1)+"""-"""+str(nPW+mPW)+""";
      partition dtype = 2:sc, cc;
      set partition = dtype;
      unlink Statefreq=(all) shape=(all);
      lset applyto=(all) rates=gamma;
      lset applyto=(1) coding=noabsencesites;
      lset applyto=(2) coding=noabsencesites;
      constraint albanian = ALBANIAN_GHEG ALBANIAN ALBANIAN_TOSK;
      constraint romance = 4-.;
      prset topologypr = constraints(albanian,romance);
      mcmcp stoprule=yes stopval = 0.01 filename = albanoRomancePW;
      mcmcp mcmcdiagn=yes diagnfreq=5000 samplefreq=10000 burninfrac=.1;
      set seed=12345;
      set swapseed=12345;
      mcmc ngen = 11000000;
end;"""

mbCommandsSW = """#NEXUS
begin MrBayes;
      execute albanoRomance_sc_cc_SW.nex;
      charset sc = 1-"""+str(nSW)+""";
      charset cc = """+str(nSW+1)+"""-"""+str(nSW+mSW)+""";
      partition dtype = 2:sc, cc;
      set partition = dtype;
      unlink Statefreq=(all) shape=(all);
      lset applyto=(all) rates=gamma;
      lset applyto=(1) coding=noabsencesites;
      lset applyto=(2) coding=noabsencesites;
      constraint albanian = ALBANIAN_GHEG ALBANIAN ALBANIAN_TOSK;
      constraint romance = 4-.;
      prset topologypr = constraints(albanian,romance);
      mcmcp stoprule=yes stopval = 0.01 filename = albanoRomanceSW;
      mcmcp mcmcdiagn=yes diagnfreq=5000 samplefreq=10000 burninfrac=.1;
      set seed=12345;
      set swapseed=12345;
      mcmc ngen = 11000000;
end;"""

with open('albanoRomance.mbPW.nex','w') as f:
    f.write(mbCommandsPW)

with open('albanoRomance.mbSW.nex','w') as f:
    f.write(mbCommandsSW)
