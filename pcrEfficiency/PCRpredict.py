import bioframe, pandas, numpy, Bio.Seq, os, random, pyranges
from amplicon import *
from predict import *

def infer_seq(ext, BioPrimerRegion, sample, path, seqlen=150, explen=300):
    # save the genome size to ref2len
    ref2len = {}
    for row in bioframe.read_table('/home/ljw/hg19_with_bowtie2_index/hg19.chrom.sizes', schema=['chrom', 'length']).itertuples():
        ref2len[row.chrom] = row.length
    genome = bioframe.load_fasta('/home/ljw/hg19_with_random/hg19.fa')
    # extract .fa or .fasta (assume only 1), path must be full so as to use os.listdir
    fastafile = [file for file in os.listdir(path) if file.endswith('.fa') or file.endswith('.fasta')][0]
    fasta = bioframe.load_fasta(os.path.join(path, fastafile))
    # extract bioprimer (strand sensitive)
    bioprimer = genome[BioPrimerRegion[0]].ff.fetch(BioPrimerRegion[0], BioPrimerRegion[1], BioPrimerRegion[2]).upper()
    if BioPrimerRegion[3] == '-':
        bioprimer = Bio.Seq.Seq(bioprimer).reverse_complement().__str__()
    queries1, queries2 = [], []
    for r in ext.itertuples():
        # extract read by name
        fa = fasta[r.QNAME].ff.fetch(r.QNAME).upper()
        # regionFirst is downstream to BioPrimerRegion
        downFirst = r.regionFirst[3]==BioPrimerRegion[3] and (BioPrimerRegion[3]=='+' and r.regionFirst[1]>=BioPrimerRegion[2] or BioPrimerRegion[3]=='-' and r.regionFirst[2]<=BioPrimerRegion[1])
        # require that the regionFirst of read (fa) must be downstream to BioPrimerRegion, and the unmapped dangle part of read must coincide the sequenced part of forward primer
        if downFirst and r.dangle<=len(sample['forwardPrimer']['seqed']) and fa[:r.dangle]==sample['forwardPrimer']['seqed'][:r.dangle]:
            # extract seq as between the end of bio-primer and the start of the first mapping region of the read (strand sensitive)
            if BioPrimerRegion[3] == '+':
                seq = genome[BioPrimerRegion[0]].ff.fetch(BioPrimerRegion[0], BioPrimerRegion[2], r.regionFirst[1]).upper()
            else:
                seq = genome[BioPrimerRegion[0]].ff.fetch(BioPrimerRegion[0], r.regionFirst[2], BioPrimerRegion[1]).upper()
                seq = Bio.Seq.Seq(seq).reverse_complement().__str__()
            # replace N in seq by random ATCG, and junc bio-primer and seq
            seq = bioprimer + randSub(seq)
            # replace N in read (fa) by random ATCG
            fa = randSub(fa)
            if len(fa) == seqlen: # fa has maximal NGS length, which implies that fa may not be fully sequenced
                # extend read from the end of the last region to the bio-primer expect length (explen, strand sensitive)
                if r.regionLast[3] == '+':
                    seg = genome[r.regionLast[0]].ff.fetch(r.regionLast[0], r.regionLast[2], min(ref2len[r.regionLast[0]],r.regionLast[2]+explen-len(fa))).upper()
                else:
                    seg = genome[r.regionLast[0]].ff.fetch(r.regionLast[0], max(0,r.regionLast[1]-(explen-len(fa))), r.regionLast[1]).upper()
                    seg = Bio.Seq.Seq(seg).reverse_complement().__str__()
                # strip downstream N
                seg = seg.rstrip('N')
                # sub random ATCG to the left Ns
                seg = randSub(seg)
            else: # fa is fully sequenced, thereby not extended
                    seg = ''

            query1 = Amplicon()
            query1.setLabel(r.QNAME)
            # PCR = BioPrimer + BioEnd2FirstStart + ReadFromMapped + Extend
            query1.setSequence(seq + fa[r.dangle:] + seg)
            # no 3' primer
            query1.setPrimerPair(bioprimer,'')
            # checking = query1.checkHybridization(query1.getSequence(), query1.getForward(),query1.getReverse())
            # query1.setSequence(checking[1])
            queries1.append(query1)

            query2 = Amplicon()
            query2.setLabel(r.QNAME)
            # remind that the dangle part of fa coinsides the sequenced part of forward primer, adapter merging reverse primer is actually the extension
            query2.setSequence(sample['forwardPrimer']['unseqed'] + fa + seg + sample['reversePrimer']['merge'])
            # the reverse primer in sample is saved as the main strand, RC it to use Amplicon
            reversePrimer = Bio.Seq.Seq(sample['reversePrimer']['primer'])
            query2.setPrimerPair(sample['forwardPrimer']['primer'], reversePrimer.reverse_complement().__str__())
            # checking = query2.checkHybridization(query2.getSequence(), query2.getForward(),query2.getReverse())
            # query2.setSequence(checking[1])
            queries2.append(query2)

    return queries1, queries2
            



def read_ext(extfile, BioPrimerRegion):
    # transform .ext file to the following dataframe:
    # QNAME \t INDEX \t regionFirst \t regionLast \t dangle
    # regionFirst is the first mapped seg for each QNAME and INDEX, which prefers downstream to BioPrimerRegion, and then prefer closer 
    # dangle is the upstream unmapped length before regionFirst
    # if query has more than one mapped segs, regionLast is that closest to BioPrimerRegion other than the first seg (not necessary downstream)
    # regionLast is regionFirst if there is only one mapped seg
    chromsizes = bioframe.read_table('/home/ljw/hg19_with_random/hg19.chrom.sizes', schema=['chrom','length'])
    qnames, indices, regionsFirst, regionsLast, dangles = [], [], [], [], []
    with open(extfile, 'r') as f:
        lines = f.readlines()
        if len(lines)==0:
            return pandas.DataFrame({'QNAME' : [], 'INDEX' : [], 'regionFirst' : [], 'regionLast' : [], 'dangle': []})
        i = 0
        while i < len(lines):
            if lines[i].startswith('>'):
                if i > 0:
                    regionsLast.append(bestRegion)
                qname, index = lines[i].split()
                qname = qname.lstrip('>')
                flagFirst = True
                minDis = numpy.inf
                FirstMinDis = numpy.inf
                i += 1
            else:
                fields1 = lines[i].split(':')[1].split()
                fields2 = lines[i+1].split(':')[1].split()
                for j in range(0, len(fields1)-2, 2):
                    if fields1[j].endswith('_RC'):
                        strand = '-'
                        chr = fields1[j][1:-3].lstrip('>')
                    else:
                        strand = '+'
                        chr = fields1[j][1:]
                    chrlen = chromsizes.loc[chromsizes['chrom']==chr,'length'].values[0]
                    if strand == '+':
                        start = int(fields1[j+1])
                        end = int(fields2[j+1])
                    elif strand == '-':
                        start = chrlen - int(fields2[j+1])
                        end = chrlen - int(fields1[j+1])
                    region = (chr, start, end, strand)
                    dis = bioframe.closest(bioframe.make_viewframe(BioPrimerRegion[:3]), bioframe.make_viewframe(region[:3]), return_input=False).values.item()
                    if flagFirst:
                        if region[3]!=BioPrimerRegion[3] or BioPrimerRegion[3]=='+' and region[1]<BioPrimerRegion[2] or BioPrimerRegion[3]=='-' and region[2]>BioPrimerRegion[1]: # region is not downstream to BioPrimerRegion
                            dis = numpy.inf
                        if dis <= FirstMinDis:
                            FirstMinDis = dis
                            bestRegion = region
                    else:
                        if dis <= minDis:
                            minDis = dis
                            bestRegion = region                    
                if flagFirst:
                    qnames.append(qname)
                    indices.append(index)
                    regionsFirst.append(bestRegion)
                    dangles.append(int(fields1[-2]))
                    flagFirst = False
                i += 2
        regionsLast.append(bestRegion)

    return pandas.DataFrame({'QNAME' : qnames, 'INDEX' : indices, 'regionFirst' : regionsFirst, 'regionLast' : regionsLast, 'dangle': dangles})

def randSub(seq):
    vc = ['A', 'T', 'C', 'G']
    seq = [*seq]
    for i in range(len(seq)):
        if seq[i] not in vc:
            seq[i] = vc[random.randint(0,3)]
    return ''.join(seq)