#!/bin/pyhton3

import pandas as pd

data_cols = ['counts', 'percentage', ]
zero_to_hundred = pd.Series(range(0,101))


# thresh_0.80:
data1 = pd.read_table('20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts', sep=' ', header=None, names=data_cols)
data2 = pd.read_table('20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts', sep=' ', header=None, names=data_cols)
data3 = pd.read_table('20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts', sep=' ', header=None, names=data_cols)
data4 = pd.read_table('20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts', sep=' ', header=None, names=data_cols)
data5 = pd.read_table('Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.80.5mC.counts', sep=' ', header=None, names=data_cols)

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
df5 = pd.DataFrame(data5)

df1.set_index("percentage", inplace=True)
df2.set_index("percentage", inplace=True)
df3.set_index("percentage", inplace=True)
df4.set_index("percentage", inplace=True)
df5.set_index("percentage", inplace=True)

df1_mod = df1.reindex(zero_to_hundred, fill_value=0)
df2_mod = df2.reindex(zero_to_hundred, fill_value=0)
df3_mod = df3.reindex(zero_to_hundred, fill_value=0)
df4_mod = df4.reindex(zero_to_hundred, fill_value=0)
df5_mod = df5.reindex(zero_to_hundred, fill_value=0)

final = pd.concat([df1_mod,df2_mod,df3_mod,df4_mod,df5_mod], axis=1, sort=False)

final.to_csv('../summaries/Aphrocallistes_megalodon_methylation_summary_thres_0.80.5mC.counts', sep="\t", index=True)


# thresh_0.75:
data1 = pd.read_table('20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts', sep=' ', header=None, names=data_cols)
data2 = pd.read_table('20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts', sep=' ', header=None, names=data_cols)
data3 = pd.read_table('20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts', sep=' ', header=None, names=data_cols)
data4 = pd.read_table('20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts', sep=' ', header=None, names=data_cols)
data5 = pd.read_table('Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.75.5mC.counts', sep=' ', header=None, names=data_cols)

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
df5 = pd.DataFrame(data5)

df1.set_index("percentage", inplace=True)
df2.set_index("percentage", inplace=True)
df3.set_index("percentage", inplace=True)
df4.set_index("percentage", inplace=True)
df5.set_index("percentage", inplace=True)

df1_mod = df1.reindex(zero_to_hundred, fill_value=0)
df2_mod = df2.reindex(zero_to_hundred, fill_value=0)
df3_mod = df3.reindex(zero_to_hundred, fill_value=0)
df4_mod = df4.reindex(zero_to_hundred, fill_value=0)
df5_mod = df5.reindex(zero_to_hundred, fill_value=0)

final = pd.concat([df1_mod,df2_mod,df3_mod,df4_mod,df5_mod], axis=1, sort=False)

final.to_csv('../summaries/Aphrocallistes_megalodon_methylation_summary_thres_0.75.5mC.counts', sep="\t", index=True)




# thresh_0.70:
data1 = pd.read_table('20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts', sep=' ', header=None, names=data_cols)
data2 = pd.read_table('20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts', sep=' ', header=None, names=data_cols)
data3 = pd.read_table('20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts', sep=' ', header=None, names=data_cols)
data4 = pd.read_table('20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts', sep=' ', header=None, names=data_cols)
data5 = pd.read_table('Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.70.5mC.counts', sep=' ', header=None, names=data_cols)

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
df5 = pd.DataFrame(data5)

df1.set_index("percentage", inplace=True)
df2.set_index("percentage", inplace=True)
df3.set_index("percentage", inplace=True)
df4.set_index("percentage", inplace=True)
df5.set_index("percentage", inplace=True)

df1_mod = df1.reindex(zero_to_hundred, fill_value=0)
df2_mod = df2.reindex(zero_to_hundred, fill_value=0)
df3_mod = df3.reindex(zero_to_hundred, fill_value=0)
df4_mod = df4.reindex(zero_to_hundred, fill_value=0)
df5_mod = df5.reindex(zero_to_hundred, fill_value=0)

final = pd.concat([df1_mod,df2_mod,df3_mod,df4_mod,df5_mod], axis=1, sort=False)

final.to_csv('../summaries/Aphrocallistes_megalodon_methylation_summary_thres_0.70.5mC.counts', sep="\t", index=True)




# thresh_0.65:
data1 = pd.read_table('20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts', sep=' ', header=None, names=data_cols)
data2 = pd.read_table('20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts', sep=' ', header=None, names=data_cols)
data3 = pd.read_table('20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts', sep=' ', header=None, names=data_cols)
data4 = pd.read_table('20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts', sep=' ', header=None, names=data_cols)
data5 = pd.read_table('Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.65.5mC.counts', sep=' ', header=None, names=data_cols)

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
df5 = pd.DataFrame(data5)

df1.set_index("percentage", inplace=True)
df2.set_index("percentage", inplace=True)
df3.set_index("percentage", inplace=True)
df4.set_index("percentage", inplace=True)
df5.set_index("percentage", inplace=True)

df1_mod = df1.reindex(zero_to_hundred, fill_value=0)
df2_mod = df2.reindex(zero_to_hundred, fill_value=0)
df3_mod = df3.reindex(zero_to_hundred, fill_value=0)
df4_mod = df4.reindex(zero_to_hundred, fill_value=0)
df5_mod = df5.reindex(zero_to_hundred, fill_value=0)

final = pd.concat([df1_mod,df2_mod,df3_mod,df4_mod,df5_mod], axis=1, sort=False)

final.to_csv('../summaries/Aphrocallistes_megalodon_methylation_summary_thres_0.65.5mC.counts', sep="\t", index=True)




# thresh_0.60:
data1 = pd.read_table('20171103_1624_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts', sep=' ', header=None, names=data_cols)
data2 = pd.read_table('20171109_1116_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts', sep=' ', header=None, names=data_cols)
data3 = pd.read_table('20180308_1413_PromethION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts', sep=' ', header=None, names=data_cols)
data4 = pd.read_table('20200221_1309_MinION_megalodon_min_modbases_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts', sep=' ', header=None, names=data_cols)
data5 = pd.read_table('Aphrocallistes_MERGED_5mC_CpG_v001_modified_bases.thresh_0.60.5mC.counts', sep=' ', header=None, names=data_cols)

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
df5 = pd.DataFrame(data5)

df1.set_index("percentage", inplace=True)
df2.set_index("percentage", inplace=True)
df3.set_index("percentage", inplace=True)
df4.set_index("percentage", inplace=True)
df5.set_index("percentage", inplace=True)

df1_mod = df1.reindex(zero_to_hundred, fill_value=0)
df2_mod = df2.reindex(zero_to_hundred, fill_value=0)
df3_mod = df3.reindex(zero_to_hundred, fill_value=0)
df4_mod = df4.reindex(zero_to_hundred, fill_value=0)
df5_mod = df5.reindex(zero_to_hundred, fill_value=0)

final = pd.concat([df1_mod,df2_mod,df3_mod,df4_mod,df5_mod], axis=1, sort=False)

final.to_csv('../summaries/Aphrocallistes_megalodon_methylation_summary_thres_0.60.5mC.counts', sep="\t", index=True)

