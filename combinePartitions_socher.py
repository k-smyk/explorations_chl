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


scPWsocher = pd.read_table('dialign+original_pipeline/albanoRomanceSC-pw.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)

scSWsocher = pd.read_table('sw+original_pipeline/albanoRomanceSC-sw.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)

scNWsocher = pd.read_table('original_results/albanoRomanceSC.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)


ccPWsocher = pd.read_table('socher_albanoRomanceCC.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)

ccSWsocher = pd.read_table('socher_albanoRomanceCC.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)

ccNWsocher = pd.read_table('socher_albanoRomanceCC.nex',
                           skiprows=7,
                           skipfooter=4, engine='python',
                           index_col=0,
                           sep='\s+', header=None)

ccPW = ccPWsocher.loc[scPWsocher.index]
ccSW = ccSWsocher.loc[scSWsocher.index]
ccNW = ccNWsocher.loc[scNWsocher.index]

scCC_PW = pd.DataFrame([list(scPWsocher[1][l] + ccPW[1][l]) for l in scPWsocher.index],
                       index = scPWsocher.index)
scCC_SW = pd.DataFrame([list(scSWsocher[1][l] + ccSW[1][l]) for l in scSWsocher.index],
                       index = scSWsocher.index)
scCC_NW = pd.DataFrame([list(scNWsocher[1][l] + ccNW[1][l]) for l in scNWsocher.index],
                       index = scNWsocher.index)

nexCharOutput(scCC_PW.values,scCC_PW.index,'socher_albanoRomance_sc_cc_PW.nex',datatype='restriction')
nexCharOutput(scCC_SW.values,scCC_SW.index,'socher_albanoRomance_sc_cc_SW.nex',datatype='restriction')
nexCharOutput(scCC_NW.values,scCC_NW.index,'socher_albanoRomance_sc_cc_NW.nex',datatype='restriction')

nPW = len(scPWsocher.values[0][0])
nSW = len(scSWsocher.values[0][0])
mPW = len(ccPW.values[0][0])
mSW = len(ccSW.values[0][0])
mNW = len(ccNW.values[0][0])
nNW = len(scNWsocher.values[0][0])

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

mbCommandNW = """#NEXUS
begin MrBayes;
      execute socher_albanoRomance_sc_cc_NW.nex;
      charset sc = 1-"""+str(nNW)+""";
      charset cc = """+str(nNW+1)+"""-"""+str(nNW+mNW)+""";
      partition dtype = 2:sc, cc;
      set partition = dtype;
      unlink Statefreq=(all) shape=(all);
      lset applyto=(all) rates=gamma;
      lset applyto=(1) coding=noabsencesites;
      lset applyto=(2) coding=noabsencesites;
      constraint albanian = ALBANIAN_GHEG ALBANIAN ALBANIAN_TOSK;
      constraint romance = 4-.;
      prset topologypr = constraints(albanian,romance);
      mcmcp stoprule=yes stopval = 0.01 filename = socher_albanoRomance;
      mcmcp mcmcdiagn=yes diagnfreq=5000 samplefreq=10000 burninfrac=.1;
      set seed=12345;
      set swapseed=12345;
      mcmc ngen = 11000000;
end;"""

# with open('dialign+original_pipeline/albanoRomance.mbPW.nex', 'w') as f:
#     f.write(mbCommandsPW)
# with open('sw+original_pipeline/albanoRomance.mbSW.nex', 'w') as f:
#     f.write(mbCommandsSW)
with open('socher_albanoRomance.mbNW.nex', 'w') as f:
    f.write(mbCommandNW)
