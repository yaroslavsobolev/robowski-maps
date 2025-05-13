
import pandas as pd
df_results = pd.read_pickle(f'hntz_df_results_model.pkl')
import pickle
pickle.HIGHEST_PROTOCOL = 4
df_results.to_hdf('hntz_df_results_model_p4.hdf', key='df')