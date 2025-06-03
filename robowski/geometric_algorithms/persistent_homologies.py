import pandas as pd
import numpy as np
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt

substances = ['ic001','am001','ald001','ptsa']
yield_column_name = 'yield'
data_file = 'postprocessed_yields_decimated.csv'

# yield_threshold = 0.11
yield_threshold = 0.07

# filter the data to only include rows where the yield is above the threshold
df_results = pd.read_csv(data_file)
df_results = df_results[df_results[yield_column_name] > yield_threshold]

# use substrate columns are coordinates
data = df_results[substances].to_numpy()

# data = np.random.random((100,2))
print(f'Data shape: {data.shape}')

diagrams = ripser(data, maxdim=2)['dgms']
plot_diagrams(diagrams, show=False)

# save the diagram figure to a file
plt.savefig(f'persistence_diagram_UGi_{yield_threshold:.5f}.png', dpi=300, bbox_inches='tight')
plt.show()