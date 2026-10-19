import rdata, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import os; D=os.path.expanduser("~/Desktop/WB Data Analysis /Phyloseq Files/")  # requires: pip install rdata
RANKS=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
def load(core):
    o=rdata.conversion.convert(rdata.parser.parse_file(D+f'{core}_phyloseq_collapsed.rds'))
    ot=o.otu_table; X=pd.DataFrame(ot.values, index=ot.coords['dim_0'].values, columns=ot.coords['dim_1'].values)
    tt=np.array(o.tax_table,dtype=object).reshape(len(RANKS),-1).T if o.tax_table.size==len(X)*7 else None
    T=pd.DataFrame(tt,index=X.index,columns=RANKS)
    return X,T
