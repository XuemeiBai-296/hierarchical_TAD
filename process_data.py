import math

import pandas as pd

stage_dict = {"8cell": 0, "ICM": 1, "mESC": 2, "late2cell": 3}
label_dict = {"raise": 0, "reduce": 1}
boundary_data_path = "boundary_data/GSE82185_17_Nature_Xiewei/"
ep_data_path = "ep_data/epigenetics/GSE82185_17_Nature_Xiewei/20kb"

cell8_ICM_boudry = pd.read_csv(
    "boundary_data/GSE82185_17_Nature_Xiewei/8cell-ICM/8cell-ICM_boundary_raise&reduce.txt"
    , sep="\t")

ICM_mESC_boudry = pd.read_csv(
    "boundary_data/GSE82185_17_Nature_Xiewei/ICM-mESC/ICM-mESC_boundary_raise&reduce.txt"
    , sep="\t")

late2cell_8cell_boudry = pd.read_csv(
    "boundary_data/GSE82185_17_Nature_Xiewei/late2cell-8cell/late2cell-8cell_boundary_raise&reduce.txt"
    , sep="\t")

cell8_ICM_atac = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/ATAC_overlap/8cell-ICM_overlap.txt"
                               , sep="\t")

ICM_mESC_atac = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/ATAC_overlap/ICM-mESC_overlap.txt"
                               , sep="\t")
late2cell_8cell_atac = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/ATAC_overlap/late2cell-8cell_overlap.txt"
                                 , sep="\t")


cell8_ICM_k4me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K4me3_overlap/8cell-ICM_overlap.txt"
                                , sep="\t")

ICM_mESC_k4me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K4me3_overlap/ICM-mESC_overlap.txt"
                                , sep="\t")

late2cell_8cell_k4me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K4me3_overlap/late2cell-8cell_overlap.txt"
                                  , sep="\t")

cell8_ICM_k9me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K9me3_overlap/8cell-ICM_overlap.txt"
                                , sep="\t")

ICM_mESC_k9me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K9me3_overlap/ICM-mESC_overlap.txt"
                                , sep="\t")

late2cell_8cell_k9me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K9me3_overlap/late2cell-8cell_overlap.txt"
                                  , sep="\t")

cell8_ICM_k27me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K27me3_overlap/8cell-ICM_overlap.txt"
                                , sep="\t")

ICM_mESC_k27me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K27me3_overlap/ICM-mESC_overlap.txt"
                                , sep="\t")

late2cell_8cell_k27me3 = pd.read_csv("epigenetics/GSE82185_17_Nature_Xiewei/20kb/K27me3_overlap/late2cell-8cell_overlap.txt"
                                  , sep="\t")

all_data = []
all_label = []

# process late2cell_8cell
for i in range(len(late2cell_8cell_boudry)):
    cur_data = late2cell_8cell_boudry.iloc[i]
    label = cur_data['group']
    label = label_dict[label]
    stages = cur_data['Cell'].split("-")
    pre_stage = stage_dict[stages[0]]
    target_stage = stage_dict[stages[1]]
    pre_level = int(cur_data['pre.level'])

    chrome = cur_data['chr']
    chrome_start = cur_data['start']
    chrome_end = cur_data['end']

    # get ATAC data
    cur_atac_data = late2cell_8cell_atac[late2cell_8cell_atac['chr.x'] == cur_data['chr']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.start'] >= cur_data['start']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.end'] <= cur_data['end']]
    atact_value = cur_atac_data['value'].mean()
    if math.isnan(atact_value):
        atact_value = 0
        print(atact_value)

    # get k4me3 data
    cur_k4me3_data = late2cell_8cell_k4me3[late2cell_8cell_k4me3['chr.x'] == cur_data['chr']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.start'] >= cur_data['start']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.end'] <= cur_data['end']]
    k4me3_value = cur_k4me3_data['value'].mean()
    if math.isnan(cur_k4me3_data['value'].mean()):
        k4me3_value = 0
        print(k4me3_value)

    # get k9me3 data
    cur_k9me3_data = late2cell_8cell_k9me3[late2cell_8cell_k9me3['chr.x'] == cur_data['chr']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.start'] >= cur_data['start']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.end'] <= cur_data['end']]
    k9me3_value = cur_k9me3_data['value'].mean()
    if math.isnan(k9me3_value):
        k9me3_value = 0
        print(k9me3_value)

    # get k27me3 data
    cur_k27me3_data = late2cell_8cell_k27me3[late2cell_8cell_k27me3['chr.x'] == cur_data['chr']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.start'] >= cur_data['start']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.end'] <= cur_data['end']]
    k27me3_value = cur_k27me3_data['value'].mean()
    if math.isnan(k27me3_value):
        k27me3_value = 0
        print(k27me3_value)

    all_data.append((pre_stage, target_stage, pre_level, atact_value, k4me3_value, k9me3_value, k27me3_value))
    all_label.append(label)

# process ICM_mESC
for i in range(len(ICM_mESC_boudry)):
    cur_data = ICM_mESC_boudry.iloc[i]
    label = cur_data['group']
    label = label_dict[label]
    stages = cur_data['Cell'].split("-")
    pre_stage = stage_dict[stages[0]]
    target_stage = stage_dict[stages[1]]
    pre_level = int(cur_data['pre.level'])

    chrome = cur_data['chr']
    chrome_start = cur_data['start']
    chrome_end = cur_data['end']

    # get ATAC data
    cur_atac_data = ICM_mESC_atac[ICM_mESC_atac['chr.x'] == cur_data['chr']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.start'] >= cur_data['start']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.end'] <= cur_data['end']]
    atact_value = cur_atac_data['value'].mean()
    if math.isnan(atact_value):
        atact_value = 0
        print(atact_value)

    # get k4me3 data
    cur_k4me3_data = ICM_mESC_k4me3[ICM_mESC_k4me3['chr.x'] == cur_data['chr']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.start'] >= cur_data['start']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.end'] <= cur_data['end']]
    k4me3_value = cur_k4me3_data['value'].mean()
    if math.isnan(k4me3_value):
        k4me3_value = 0
        print(k4me3_value)

    # get k9me3 data
    cur_k9me3_data = ICM_mESC_k9me3[ICM_mESC_k9me3['chr.x'] == cur_data['chr']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.start'] >= cur_data['start']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.end'] <= cur_data['end']]
    k9me3_value = cur_k9me3_data['value'].mean()
    if math.isnan(k9me3_value):
        k9me3_value = 0
        print(k9me3_value)

    # get k27me3 data
    cur_k27me3_data = ICM_mESC_k27me3[ICM_mESC_k27me3['chr.x'] == cur_data['chr']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.start'] >= cur_data['start']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.end'] <= cur_data['end']]
    k27me3_value = cur_k27me3_data['value'].mean()
    if math.isnan(k27me3_value):
        k27me3_value = 0
        print(k27me3_value)

    all_data.append((pre_stage, target_stage, pre_level, atact_value, k4me3_value, k9me3_value, k27me3_value))
    all_label.append(label)

# process late2cell_8cell
for i in range(len(late2cell_8cell_boudry)):
    cur_data = late2cell_8cell_boudry.iloc[i]
    label = cur_data['group']
    label = label_dict[label]
    stages = cur_data['Cell'].split("-")
    pre_stage = stage_dict[stages[0]]
    target_stage = stage_dict[stages[1]]
    pre_level = int(cur_data['pre.level'])

    chrome = cur_data['chr']
    chrome_start = cur_data['start']
    chrome_end = cur_data['end']

    # get ATAC data
    cur_atac_data = late2cell_8cell_atac[late2cell_8cell_atac['chr.x'] == cur_data['chr']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.start'] >= cur_data['start']]
    cur_atac_data = cur_atac_data[cur_atac_data['TAD.end'] <= cur_data['end']]
    atact_value = cur_atac_data['value'].mean()
    if math.isnan(atact_value):
        atact_value = 0
        print(atact_value)

    # get k4me3 data
    cur_k4me3_data = late2cell_8cell_k4me3[late2cell_8cell_k4me3['chr.x'] == cur_data['chr']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.start'] >= cur_data['start']]
    cur_k4me3_data = cur_k4me3_data[cur_k4me3_data['TAD.end'] <= cur_data['end']]
    k4me3_value = cur_k4me3_data['value'].mean()
    if math.isnan(k4me3_value):
        k4me3_value = 0
        print(k4me3_value)

    # get k9me3 data
    cur_k9me3_data = late2cell_8cell_k9me3[late2cell_8cell_k9me3['chr.x'] == cur_data['chr']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.start'] >= cur_data['start']]
    cur_k9me3_data = cur_k9me3_data[cur_k9me3_data['TAD.end'] <= cur_data['end']]
    k9me3_value = cur_k9me3_data['value'].mean()
    if math.isnan(k9me3_value):
        k9me3_value = 0
        print(k9me3_value)

    # get k27me3 data
    cur_k27me3_data = late2cell_8cell_k27me3[late2cell_8cell_k27me3['chr.x'] == cur_data['chr']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.start'] >= cur_data['start']]
    cur_k27me3_data = cur_k27me3_data[cur_k27me3_data['TAD.end'] <= cur_data['end']]
    k27me3_value = cur_k27me3_data['value'].mean()
    if math.isnan(k27me3_value):
        k27me3_value = 0
        print(k27me3_value)

    all_data.append((pre_stage, target_stage, pre_level, atact_value, k4me3_value, k9me3_value, k27me3_value))
    all_label.append(label)

import pickle

pickle.dump(all_data, open('all_data.pkl', 'wb'))
pickle.dump(all_label, open('all_label.pkl', 'wb'))
