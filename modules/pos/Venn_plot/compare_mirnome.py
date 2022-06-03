#!/usr/bin/env python3

""" This scrip take the TARGET_IDs data and make venn Diagram

Basic ideas:
    1. Read list with the IDs
    2. Plot venn Diagram

Autor:  Eduardo García López"""

## Import python libraries
import sys
import subprocess
import pandas as  pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3

## Read args from command line
    ## Uncomment For debugging only
    ## Comment for production mode only
#sys.argv = ("0", "All_targets.ref.tsv", "All_targets.alt.tsv")

##get targets files
targets_ref = sys.argv[1]

targets_alt = sys.argv[2]

#Read file as list
ref_list = pd.read_csv(targets_ref, sep = '\t')
alt_list = pd.read_csv(targets_alt, sep = '\t')


ref_targets = ref_list.groupby('prediction_tool').get_group('both')

alt_targets = alt_list.groupby('prediction_tool').get_group('both')


venn2([set(ref_targets['target_ID']),
       set(alt_targets['target_ID'])],
       set_labels=('Ref_targets', 'Alt_targets')
     )

plt.title('Venn Diagram')

#plt.show()

plt.savefig('miRNome_changes.png', dpi=(800))

plt.close()
