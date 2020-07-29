'''
Created on Feb 19, 2014

@author: yonatanf
'''
from numpy import genfromtxt, all, abs


def read_SparCC_output(path):
    '''
    Read one of SparCC's output files (correlations or p-values), and return a tuple of labels and an array of values.
    '''
    tmp    = genfromtxt(path, names=True, excludelist=['OTU_id'])
    labels = tmp.dtype.names[1:]
    vals   = tmp.view((float, len(tmp.dtype.names)))[:,1:]
    return vals, labels

def main(cor_path, pval_path, out_path, cor_th, pval_th, delim='\t'):
    ## read data
    cors, otus   = read_SparCC_output(cor_path)
    pvals, potus = read_SparCC_output(pval_path)
    if not all(otus==potus):
        raise ValueError, 'OTU labels in correlation and p-value files do not match!'
    ## go over all pairs and find significant ones
    n = len(otus)
    lines =  delim.join(['OTU 1', 'OTU 2', 'corr', 'pval']) +'\n'
    for i in xrange(n):
        for j in xrange(i+1,n):
            c = cors[i,j]
            p = pvals[i,j]
            if abs(c)>=cor_th and p<=pval_th: ## Modify this line to change the filtering criteria
                o1 = otus[i]
                o2 = otus[j]
                line = delim.join([o1, o2, str(c), str(p)]) +'\n'
                lines += line
    # write lines to output file
    with open(out_path, 'w') as fw:
        fw.write(lines)
        

def test_read_SparCC_output():
    path = 'pvals_one_sided.txt'
    vals, labels = read_SparCC_output(path)
    print labels[:3]
    print vals[:3,:3]
    print len(labels), vals.shape

if __name__ == '__main__':
    cor_path  = 'cor_sparcc.out'
    pval_path = 'pvals_two_sided.txt'
    out_path  = 'selected_cor.txt'
    cor_th    = 0.6
    pval_th   = 0.05
    main(cor_path, pval_path, out_path, cor_th, pval_th)

