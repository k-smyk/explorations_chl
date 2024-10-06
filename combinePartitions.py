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


scPW = pd.read_table('dialign+original_pipeline/albanoRomanceSC-pw.nex',
					 skiprows=7,
					 skipfooter=4, engine='python',
					 index_col=0,
					 sep='\s+', header=None)

scSW = pd.read_table('sw+original_pipeline/albanoRomanceSC-sw.nex',
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

ccPWsocher = ccPWsocher.loc[scPW.index]
ccSWsocher = ccSWsocher.loc[scSW.index]

scCC_PW = pd.DataFrame([list(scPW[1][l] + ccPWsocher[1][l]) for l in scPW.index],
					   index = scPW.index)
scCC_SW = pd.DataFrame([list(scSW[1][l] + ccSWsocher[1][l]) for l in scSW.index],
					   index = scSW.index)

nexCharOutput(scCC_PW.values,scCC_PW.index,'socher_albanoRomance_sc_cc_PW.nex',datatype='restriction')
nexCharOutput(scCC_SW.values,scCC_SW.index,'socher_albanoRomance_sc_cc_SW.nex',datatype='restriction')

nPW = len(scPW.values[0][0])
nSW = len(scSW.values[0][0])
mPW = len(ccPWsocher.values[0][0])
mSW = len(ccSWsocher.values[0][0])

mbCommandsPW = """#NEXUS
begin MrBayes;
      execute socher_albanoRomance_sc_cc_PW.nex;
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
      execute socher_lbanoRomance_sc_cc_SW.nex;
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

with open('dialign+socher+original/socher_albanoRomance.mbPW.nex', 'w') as f:
    f.write(mbCommandsPW)

with open('sw+socher+original/socher_albanoRomance.mbSW.nex', 'w') as f:
    f.write(mbCommandsSW)
