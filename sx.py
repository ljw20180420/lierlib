import sys
sys.path.append('./pcrEfficiency')
import pcrEfficiency.PCRpredict, seaborn, pandas, numpy, os, pandas, bioframe, Bio.Seq, matplotlib.pyplot, seaborn, collections

### LAMT

BioPrimerRegion = ('chrX', 133607339, 133607362, '+')
Nested_HPRT1_E2F_a = {'primer' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTACCATGCTGAGGATTTGGAAAGGG', 'seqed' : 'GCTACCATGCTGAGGATTTGGAAAGGG', 'unseqed' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'}
P7_A = {'primer' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTACGAATCTCGTATGCCGTCTTCTGCTTG', 'adapter' : 'CAGATCGGAAGAGCACACGTC', 'merge' : 'CAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTACGAATCTCGTATGCCGTCTTCTGCTTG'}
Nested_HPRT1_E2F_b = {'primer' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTGCACTGCTGAGGATTTGGAAAGGGT', 'seqed' : 'GCTGCACTGCTGAGGATTTGGAAAGGGT', 'unseqed' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'}
P7_B = {'primer' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCTGAATCTCGTATGCCGTCTTCTGCTTG', 'adapter' : 'CAGATCGGAAGAGCACACGTC', 'merge' : 'CAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCTGAATCTCGTATGCCGTCTTCTGCTTG'}
Nested_HPRT1_E2F_c = {'primer' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCTATCGTGCTGAGGATTTGGAAAGGGT', 'seqed' : 'GCTATCGTGCTGAGGATTTGGAAAGGGT', 'unseqed' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'}
P7_C = {'primer' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGCGAATCTCGTATGCCGTCTTCTGCTTG', 'adapter' : 'CAGATCGGAAGAGCACACGTC', 'merge' : 'CAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGCGAATCTCGTATGCCGTCTTCTGCTTG'}
Nested_HPRT1_E2F_d = {'primer' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGATTCATTGCTGAGGATTTGGAAAGGGT', 'seqed' : 'GATTCATTGCTGAGGATTTGGAAAGGGT', 'unseqed' : 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGATT'}
P7_D = {'primer' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTCTCGATCTCGTATGCCGTCTTCTGCTTG', 'adapter' : 'CAGATCGGAAGAGCACACGTC', 'merge' : 'CAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTCTCGATCTCGTATGCCGTCTTCTGCTTG'}
Ct1_Rep1 = Ct1_Rep2 = Kappa_Rep1 = {'forwardPrimer' : Nested_HPRT1_E2F_a, 'reversePrimer' : P7_A}
Lambda_Rep1 = Lambda_Rep2 = Delta_Rep1 = {'forwardPrimer' : Nested_HPRT1_E2F_b, 'reversePrimer' : P7_B}
Kappa_Rep2 = Mu_Rep1 = Mu_Rep2 = {'forwardPrimer' : Nested_HPRT1_E2F_c, 'reversePrimer' : P7_C}
Delta_Rep2 = Theta_Rep1 = Theta_Rep2 = {'forwardPrimer' : Nested_HPRT1_E2F_d, 'reversePrimer' : P7_D}

myPredict = pcrEfficiency.PCRpredict.predict()
names = ['Ct1-Rep1', 'Ct1-Rep2', 'Kappa-Rep1', 'Kappa-Rep2', 'Lambda-Rep1', 'Lambda-Rep2', 'Delta-Rep1', 'Delta-Rep2', 'Mu-Rep1', 'Mu-Rep2', 'Theta-Rep1', 'Theta-Rep2']
output = '/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/pcr_efficiency_results'
# names = ['Ct1-Rep1']
for name in names:
    sample = eval(name.replace('-','_'))
    path = f'/home/ljw/new_fold/old_desktop/wuqiang/shoujia/finalresults/single/{name}'
    files = [file for file in os.listdir(path) if file.endswith('.ext')]
    ext = pandas.concat([pcrEfficiency.PCRpredict.read_ext(os.path.join(path,file), BioPrimerRegion) for file in files]).convert_dtypes(convert_integer=True)
    queries1, queries2 = pcrEfficiency.PCRpredict.infer_seq(ext, BioPrimerRegion, sample, path, seqlen=150, explen=300)

    for effname, queries in zip(['eff1', 'eff2'], [queries1, queries2]):
        myPredict.writeHandleForR(queries, len(queries))
        effs = numpy.array([float(eff) for eff in myPredict.predictGam()])
        pcrEfficiency.PCRpredict.call(["rm", "gamResult.data"])
        pcrEfficiency.PCRpredict.call(["rm", "primerDataForR.dat"])

        queriesLengths = numpy.array([float(len(query.getSequence())) for query in queries])
        z_scores_effs = (effs - numpy.mean(effs)) / numpy.std(effs)
        z_scores_lengths = (queriesLengths - numpy.mean(queriesLengths)) / numpy.std(queriesLengths)
        outliers = numpy.logical_or(numpy.abs(z_scores_effs) > 5, numpy.abs(z_scores_lengths) > 5)
        
        f = seaborn.displot(x=queriesLengths[~outliers], kind="hist", kde=False, bins=30)
        f.ax.set(xlabel='sequence length')
        f.savefig(os.path.join(output, 'length', f'{name}.{effname}.length.png'))
        
        f = seaborn.displot(x=effs[~outliers], kind="hist", kde=False, bins=30)
        f.ax.set(xlabel='PCR efficiency')
        f.savefig(os.path.join(output, 'effi', f'{name}.{effname}.png'))

        effsres = effs[~outliers].copy()
        effsres[effsres>2.0] = 2.0
        effsres[effsres<1.5] = 1.5
        f = seaborn.displot(x=effsres, kind="hist", kde=False, binrange=[1,2], bins=30)
        f.ax.set(xlabel='PCR efficiency')
        f.savefig(os.path.join(output, 'restri', f'{name}.{effname}.restri.png'))

### Other
sgpr = bioframe.read_table("sgpr.bed", schema=["name", "shortname", "chrom", "sg1", "sg2", "F1u", "F1s", "R1u", "R1s", "F2u", "F2s", "R2u", "R2s"])

genome_file = '/home/ljw/hg19_with_bowtie2_index/hg19.fa'
genome = bioframe.load_fasta(genome_file)
names, chrom1s, start1s, end1s, strand1s, cut1s, rc1s, pr1us, chrom2s, start2s, end2s, strand2s, cut2s, rc2s, pr2us = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
sg1skip, sg2skip = 4, 4
for row in sgpr.itertuples():
    chromseq = genome[row.chrom].ff.fetch(row.chrom).upper()
    names.extend([row.name+"_"+sfx for sfx in ["DEL","VF","VD","DUP"]])
    chrom1s.extend([row.chrom] * 4)
    chrom2s.extend([row.chrom] * 4)
    if row.name=="PRDM5":
        sg1 = Bio.Seq.Seq(row.sg1[sg1skip:]).reverse_complement().__str__()
        cut1s.extend([chromseq.find(sg1) + 3] * 4)
    else:
        sg1 = row.sg1[sg1skip:]
        cut1s.extend([chromseq.find(sg1) + len(sg1) - 3] * 4)
    if row.name=="YY1":
        sg2 = row.sg2[sg2skip:]
        cut2s.extend([chromseq.find(sg2) + len(sg2) - 3] * 4)
    else:
        sg2 = Bio.Seq.Seq(row.sg2[sg2skip:]).reverse_complement().__str__()
        cut2s.extend([chromseq.find(sg2) + 3] * 4)
    F1start, R1start, F2start, R2start = chromseq[:cut1s[-1]].rfind(row.F1s), cut1s[-1]+chromseq[cut1s[-1]:].find(Bio.Seq.Seq(row.R1s).reverse_complement().__str__()), chromseq[:cut2s[-1]].rfind(row.F2s), cut2s[-1]+chromseq[cut2s[-1]:].find(Bio.Seq.Seq(row.R2s).reverse_complement().__str__())
    start1s.extend([F1start,F1start,R1start,R1start])
    start2s.extend([R2start,F2start,R2start,F2start])
    end1s.extend([F1start+len(row.F1s)]*2+[R1start+len(row.R1s)]*2)
    end2s.extend([R2start+len(row.R2s), F2start+len(row.F2s)]*2)
    strand1s.extend(["+", "+", "-", "-"])
    strand2s.extend(["-", "+", "-", "+"])
    rc1s.extend([False, False, True, False])
    rc2s.extend([False, True, False, False])
    pr1us.extend([row.F1u, row.F1u, row.R1u, row.R1u])
    pr2us.extend([row.R2u, row.F2u, row.R2u, row.F2u])

primerPairs = pandas.DataFrame({"name" : names, "chrom1" : chrom1s, "start1" : start1s, "end1" : end1s, "strand1" : strand1s, "cut1" : cut1s, "rc1" : rc1s, "pr1u" : pr1us, "chrom2" : chrom2s, "start2" : start2s, "end2" : end2s, "strand2" : strand2s, "cut2" : cut2s, "rc2" : rc2s, "pr2u" : pr2us})

queries = []
for row in primerPairs.itertuples():
    chromseq = genome[row.chrom1].ff.fetch(row.chrom2).upper()
    print(row.name)
    if row.strand1=="+" and row.cut1<row.end1 or row.strand1=="-" and row.cut1>row.start1:
        raise Exception("primer1 is cut")
    if row.strand2=="+" and row.cut2<row.end2 or row.strand2=="-" and row.cut2>row.start2:
        raise Exception("primer2 is cut")
    strandrc1 = row.strand1 if not row.rc1 else "+-".replace(row.strand1,"")
    strandrc2 = row.strand2 if not row.rc2 else "+-".replace(row.strand2,"")
    if strandrc1==strandrc2:
        raise Exception("two primers are on the same strand")
    primers, seqs = [], []
    for start, end, strand, cut, rc, pru in zip([row.start1,row.start2],[row.end1,row.end2],[row.strand1,row.strand2],[row.cut1,row.cut2],[row.rc1,row.rc2],[row.pr1u,row.pr2u]):
        primers.append(chromseq[start:end])
        seqs.append(chromseq[min(start,end,cut):max(start,end,cut)])
        if strand=="-":
            primers[-1] = Bio.Seq.Seq(primers[-1]).reverse_complement().__str__()
            seqs[-1] = Bio.Seq.Seq(seqs[-1]).reverse_complement().__str__()
        primers[-1] = pru + primers[-1]
        seqs[-1] = pru + seqs[-1]

    seqs[-1] = Bio.Seq.Seq(seqs[-1]).reverse_complement().__str__()

    queries.append(pcrEfficiency.PCRpredict.Amplicon())
    queries[-1].setLabel(row.name)
    queries[-1].setSequence(seqs[0]+seqs[1])
    queries[-1].setPrimerPair(primers[0], primers[1])

myPredict = pcrEfficiency.PCRpredict.predict()
myPredict.writeHandleForR(queries, len(queries))
primerPairs["effs"] = numpy.array([float(eff) for eff in myPredict.predictGam()])
pcrEfficiency.PCRpredict.call(["rm", "gamResult.data"])
pcrEfficiency.PCRpredict.call(["rm", "primerDataForR.dat"])

f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
seaborn.barplot(data=primerPairs, x="name", y="effs", ax=ax)
for patch in ax.patches:
    ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{patch.get_height():.2f}', ha="center", va="bottom")
ax.tick_params(axis='x', rotation=90)
f.savefig("effs.rearr.pdf")
matplotlib.pyplot.close(f)

file = "K70-MA-DEL_merge.fastq.gz"
locus = sgpr[[True if file.find(shortname)!=-1 else False for shortname in sgpr["shortname"]]].reset_index(drop=True)
chromseq = genome[locus.loc[0,"chrom"]].ff.fetch(locus.loc[0,"chrom"])
if locus.loc[0,"name"]=="YY1":
    if (file.find("DEL")!=-1 or file.find("VF")!=-1):
        pr1s = locus.loc[0,"F1s"]
        pr1u = locus.loc[0,"F1u"]
    else:
        pr1s = locus.loc[0,"R1s"]
        pr1u = locus.loc[0,"R1u"]
else:
    if (file.find("DUP")!=-1 or file.find("VF")!=-1):
        pr1s = locus.loc[0,"F2s"]
        pr1u = locus.loc[0,"F2u"]
    else:
        pr1s = locus.loc[0,"R2s"]
        pr1u = locus.loc[0,"R2u"]
if locus.loc[0,"name"]=="YY1":
    if (file.find("DEL")!=-1 or file.find("VD")!=-1):
        pr2s = locus.loc[0,"R2s"]
        pr2u = locus.loc[0,"R2u"]
    else:
        pr2s = locus.loc[0,"F2s"]
        pr2u = locus.loc[0,"F2u"]
else:
    if (file.find("DEL")!=-1 or file.find("VF")!=-1):
        pr2s = locus.loc[0,"F1s"]
        pr2u = locus.loc[0,"F1u"]
    else:
        pr2s = locus.loc[0,"R1s"]
        pr2u = locus.loc[0,"R1u"]
pr2s = Bio.Seq.Seq(pr2s).reverse_complement().__str__()
pr2u = Bio.Seq.Seq(pr2u).reverse_complement().__str__()
f = open("/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/293T-0531/re1/merge/K70-MA-DEL_merge.fastq.gz", 'r')
duplines = collections.Counter(f.readlines()[1::4])
ids, queries = [], []
for i, seq_c in enumerate(duplines.most_common(), start=1):
    dec_seq = seq_c[0].decode()[:-1]
    pr1spos, pr2spos = dec_seq.find(pr1s), dec_seq.rfind(pr2s)
    if pr1spos!=-1 and pr2spos!=-1:
        ids.append(i)
        queries.append(pcrEfficiency.PCRpredict.Amplicon())
        queries[-1].setLabel(f"{i}")
        queries[-1].setSequence(pr1u + dec_seq[pr1spos:pr2spos+len(pr2s)] + pr2u)
        queries[-1].setPrimerPair(pr1u+pr1s, Bio.Seq.Seq(pr2s+pr2u).reverse_complement().__str__())
f.close()