import sys
sys.path.append('./pcrEfficiency')
import pcrEfficiency.PCRpredict, seaborn, pandas, numpy, os, bioframe, Bio.Seq, matplotlib.pyplot, collections, subprocess, re

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
    if not row.name.startswith("YY1"):
        primers = primers[::-1]
        seqs = seqs[::-1]

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

primerPairs["effsnorm1"] = primerPairs["effs"] / numpy.max(primerPairs["effs"])

primerPairs.to_csv("effs.rearr.tsv", sep = "\t", index = False)

f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
seaborn.barplot(data=primerPairs, x="name", y="effs", ax=ax)
for patch in ax.patches:
    ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{patch.get_height():.2f}', ha="center", va="bottom")
ax.tick_params(axis='x', rotation=90)
f.savefig("effs.rearr.pdf")
matplotlib.pyplot.close(f)

f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
seaborn.barplot(data=primerPairs, x="name", y="effsnorm1", ax=ax)
for patch in ax.patches:
    ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{patch.get_height():.2f}', ha="center", va="bottom")
ax.tick_params(axis='x', rotation=90)
f.savefig("effs.rearr.norm1.pdf")
matplotlib.pyplot.close(f)

paths = ["/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/293T-0531/re1/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/293T-0531/re2/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/HEC-1-B-0616/re1/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/HEC-1-B-0616/re2/merge/"]

with open("std_count", 'w') as fw:
    fw.write(f"name\tmean\tvariance\tcount\n")
    for path in paths:
        for file in [fi for fi in os.listdir(path) if fi.endswith(".fastq.gz")]:
            if os.stat(os.path.join(path, file)).st_size == 0:
                # print(f"{os.path.join(path, file)} is empty")
                # fw.write(f"{os.path.join(path, file)}\t{numpy.nan}\t{numpy.nan}\t0\n")
                continue
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
            f = open(os.path.join(path,file), 'r')
            duplines = collections.Counter(f.readlines()[1::4])
            ids, queries, counts = [], [], []
            for i, seq_c in enumerate(duplines.most_common()):
                counts.append(seq_c[1])
                pr1spos, pr2spos = seq_c[0].find(pr1s), seq_c[0].rfind(pr2s)
                if pr1spos!=-1 and pr2spos!=-1:
                    ids.append(i)
                    queries.append(pcrEfficiency.PCRpredict.Amplicon())
                    queries[-1].setLabel(f"{i}")
                    queries[-1].setSequence(pr1u + seq_c[0][pr1spos:pr2spos+len(pr2s)] + pr2u)
                    queries[-1].setPrimerPair(pr1u+pr1s, Bio.Seq.Seq(pr2s+pr2u).reverse_complement().__str__())

            df = pandas.DataFrame({"count" : counts})
            df["eff"] = numpy.nan
            myPredict = pcrEfficiency.PCRpredict.predict()
            myPredict.writeHandleForR(queries, len(queries))
            df.loc[ids,"eff"] = numpy.array([float(eff) for eff in myPredict.predictGam()])
            pcrEfficiency.PCRpredict.call(["rm", "gamResult.data"])
            pcrEfficiency.PCRpredict.call(["rm", "primerDataForR.dat"])
            df.to_csv(f"{os.path.join(path,file)}.eff.count", sep='\t', header=True, index=False)

            mean_val = numpy.average(df["eff"][~df["eff"].isna()], weights=df["count"][~df["eff"].isna()])
            var_val = numpy.average((df["eff"][~df["eff"].isna()] - mean_val)**2, weights=df["count"][~df["eff"].isna()])
            fw.write(f'{os.path.join(path, file)}\t{mean_val}\t{var_val}\t{numpy.sum(df["count"][~df["eff"].isna()])}\n')

            f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
            seaborn.histplot(data=df, x="eff", weights="count", binwidth=0.01, stat="percent", ax=ax)
            for patch in ax.patches:
                ax.text(patch.get_x() + patch.get_width()/2., patch.get_height(), f'{patch.get_height():.5f}', rotation=90, ha="center", va="bottom")
            ax.tick_params(axis='x', rotation=90)
            f.savefig(f"{os.path.join(path,file)}.effs.hist.pdf")
            matplotlib.pyplot.close(f)


# paths = ["/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/293T-0531/re1/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/293T-0531/re2/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/HEC-1-B-0616/re1/merge/", "/media/ljw/029ca15e-f998-413c-984e-3b63921e5f6a1/HEC-1-B-0616/re2/merge/"]
# dfs = []
# for path in paths:
#     for file in os.listdir(path):
#         if not file.endswith(".fastq.gz") or os.stat(os.path.join(path, file)).st_size == 0:
#             continue
#         df = pandas.read_csv(f"{os.path.join(path,file)}.eff.count", sep='\t')
#         df["name"] = os.path.join(path, file)
#         dfs.append(df)
# dfall = pandas.concat(dfs)
# f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
# seaborn.barplot(data=dfall, x="name", y="eff", ax=ax)
# ax.tick_params(axis='x', rotation=90)
# f.savefig(f"all_sample.boxplot.pdf")
# matplotlib.pyplot.close(f)

df_summary = pandas.read_csv("std_count", sep='\t')
loc_re = []
for name in df_summary["name"]:
    loc_re += re.findall(r"-(?:ME|PR|YY|PA|MA)(?:-)(?:DEL|VF|VR|DUP)_", name)
df_summary["loc_re"] = loc_re
df_summary["norm"] = df_summary["mean"] / max(df_summary["mean"])
data = df_summary.sort_values(by=["loc_re", "norm"]).reset_index(drop=True)
f, ax = matplotlib.pyplot.subplots(figsize=(30,30))
ax.bar(range(len(data)), data["norm"], color='b', alpha=0.5)
step = len(data)//20
ax.set(xticks=numpy.arange(step//2,len(data),step), xticklabels=data["loc_re"].unique())
f.savefig(f"all_sample.barplot.pdf")
matplotlib.pyplot.close(f)


### get examples

subprocess.check_output("> examples", shell=True)
path = "/home/ljw/wuqiang/sxdata"
num = 10
tofinds = ["-DCK-U-", "-DCK-D-", "-HPRT-U-", "-HPRT-D-"]
genome = bioframe.load_fasta("/home/ljw/hg19_with_bowtie2_index/hg19.fa")
for tofind in tofinds:
    subprocess.check_output(f"echo {tofind} >> examples", shell=True)
    dfs = []
    for file in os.listdir(path):
        if not file.endswith(".fastq.gz.fa") or file.find(tofind)==-1:
            continue
        df = pandas.read_csv(os.path.join(path, f"{file}.indel.bak"), sep='\t')
        df["exp"] = file
        dfs.append(df)
    df = pandas.concat(dfs).reset_index(drop=True)

    chromseq = genome[df.loc[0,'chrom1']].ff.fetch(df.loc[0,'chrom1'])

    strand = df.loc[0, "strand1"]
    for name, mask in zip(["single","delete","inverse"], [(df["cut1"]==df["cut2"]) & (df["key2"]!="*"), (df["cut1"]!=df["cut2"]) & (df["strand1"]==df["strand2"]), (df["cut1"]!=df["cut2"]) & (df["strand1"]!=df["strand2"])]):
        df_select = df[mask].reset_index(drop=True)
        df_select["key2"] = df_select["key2"].astype(int)
        df_select["mapstart2"] = df_select["mapstart2"].astype(int)
        df_select["mapend2"] = df_select["mapend2"].astype(int)
        if name=="inverse":
            mask2 = ((df_select["strand1"]=="+") & (df_select["mapstart2"]>=df_select["cut1"])) | ((df_select["strand1"]=="-") & (df_select["mapend2"]<=df_select["cut1"]))
            df_select = df_select[mask2].reset_index(drop=True)

        if strand=="+" and name in ["single","delete"] or strand=="-" and name=="inverse":
            df_select["resc"] = df_select["key2"]-df_select["cut2"]
        else:
            df_select["resc"] = df_select["cut2"]-df_select["key2"]
        if strand=="+":
            df_select["resc1"] = df_select["cut1"]-df_select["key1"]
        else:
            df_select["resc1"] = df_select["key1"]-df_select["cut1"]
        df_select = df_select[(df_select["resc"]>100) * (df_select["resc1"]>=0)]
        # df_select = df_select.sort_values(by=["resc"])[-min(len(df_select),num):]
        df_select = df_select.sample(min(len(df_select),num))

        for example in df_select.itertuples():
            ref1 = example.ref1_seg
            query1 = example.query1_seg
            if example.strand1=="+":
                ref1 += chromseq[example.key1:example.cut1+6]
                query1 += " "*(example.cut1-example.key1+6)
            else:
                ref1 += Bio.Seq.Seq(chromseq[example.cut1-6:example.key1]).reverse_complement().__str__()
                query1 += " "*(example.key1-example.cut1+6)
            if example.strand1==example.strand2:
                ref2 = example.ref2_seg
                query2 = example.query2_seg
                if example.strand2=="+":
                    ref2 = chromseq[example.cut2-6:example.key2] + ref2
                    query2 = " "*(example.key2-example.cut2+6) + query2
                else:
                    ref2 = Bio.Seq.Seq(chromseq[example.key2:example.cut2+6]).reverse_complement().__str__() + ref2
                    query2 = " "*(example.cut2-example.key2+6) + query2
                
                refprint = ref1[:6] + f"_{len(ref1)-41}_" + ref1[-35:] + f"_{example.cut2-example.cut1-12}_" + ref2[:23] + f"_{abs(example.key2-example.cut2)-17}_" + ref2[abs(example.key2-example.cut2)+6:abs(example.key2-example.cut2)+26]
                queryprint = query1[:6] + f"_{len(query1)-41}_" + query1[-35:] + f"_{example.cut2-example.cut1-12}_" + query2[:23] + f"_{abs(example.key2-example.cut2)-17}_" + query2[abs(example.key2-example.cut2)+6:abs(example.key2-example.cut2)+26]
            else:
                if example.strand1=="+":
                    ref2 = chromseq[example.cut1+6:example.mapend2]
                    query2 = " "*(example.mapstart2-example.cut1-6) + example.query2_seg[::-1]
                    len23 = example.cut2-6-example.mapend2
                    ref3 = chromseq[example.cut2-6:example.cut2+40]
                    query3 = " "*46
                else:
                    ref2 = Bio.Seq.Seq(chromseq[example.mapstart2:example.cut1-6]).reverse_complement().__str__()
                    query2 = " "*(example.cut1-6-example.mapend2) + example.query2_seg[::-1]
                    len23 = example.mapstart2 - example.cut2 - 6
                    ref3 = Bio.Seq.Seq(chromseq[example.cut2-40:example.cut2+6]).reverse_complement().__str__()
                    query3 = " "*46

                refprint = ref1[:6] + f"_{len(ref1)-41}_" + ref1[-35:] + f"_{len(ref2)-20}_" + ref2[-20:] + f"_{len23}_" + ref3
                queryprint = query1[:6] + f"_{len(query1)-41}_" + query1[-35:] + f"_{len(query2)-20}_" + query2[-20:] + f"_{len23}_" + query3

            alg = os.path.join(path,f"{example.exp}.alg")
            subprocess.check_output(f"echo {example.exp} >> examples; echo {name} >> examples", shell=True)
            subprocess.check_output(f"sed -n '/{example.qname}\t/,+1p' {os.path.join(path,example.exp)} >> examples", shell=True)
            subprocess.check_output(f"echo mapstart1:{example.mapstart1}\tmapend1:{example.mapend1}\tmapstart2{example.mapstart2}\tmapend2{example.mapend2} >> examples", shell=True)
            subprocess.check_output(f"echo {example.resc1}\t{example.resc} >> examples", shell=True)
            subprocess.check_output(f"echo '{refprint}' >> examples", shell=True)
            subprocess.check_output(f"echo '{queryprint}' >> examples", shell=True)
