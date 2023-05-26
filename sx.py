import sys
sys.path.append('./pcrEfficiency')
import pcrEfficiency.PCRpredict, seaborn, pandas, numpy, os


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