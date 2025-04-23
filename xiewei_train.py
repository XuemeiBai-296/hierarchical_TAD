import math
import pickle
from random import shuffle
import torch
from torch import nn
import torchmetrics
from model import boundry_model
from sklearn.model_selection import KFold
all_data = pickle.load(open('Xiewei_all_data.pkl', 'rb'))
all_label = pickle.load(open('Xiewei_all_label.pkl', 'rb'))

all_index = [i for i in range(len(all_data))]


# shuffle(all_index)
# train_index = all_index[:int(0.7 * len(all_index))]
# test_index = all_index[int(0.7 * len(all_index)):]
#
# train_data = [all_data[i] for i in train_index]
# test_data = [all_data[i] for i in test_index]
#
# train_label = [all_label[i] for i in train_index]
# test_label = [all_label[i] for i in test_index]


import torch.utils.data as data
from torch.utils.data import DataLoader

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, f1_score, accuracy_score
class Data(data.Dataset):
    def __init__(self, datasets, labels):
        self.datasets = datasets
        self.labels = labels

    def __len__(self):
        return len(self.datasets)

    def __getitem__(self, idx):
        pre_stage, target_stage, pre_level, atact_value, k4me3_value, k9me3_value, k27me3_value = self.datasets[idx]
        if math.isnan(atact_value):
            print("atac")
            print(idx)
            exit(0)
        if math.isnan(k4me3_value):
            print("k4me3_value")
            print(idx)
            exit(0)
        if math.isnan(k9me3_value):
            print("k9me3_value")
            print(idx)
            exit(0)
        if math.isnan(k27me3_value):
            print("k27me3_value")
            print(idx)
            exit(0)
        label = self.labels[idx]
        return pre_stage, target_stage, pre_level, atact_value, k4me3_value, k9me3_value, k27me3_value, label

batch_size = 256

acc_metric = torchmetrics.classification.Accuracy(task="multiclass", num_classes=2)
auc_metric = torchmetrics.classification.AUROC(task="multiclass", num_classes=2)

kf = KFold(n_splits=10, shuffle=True, random_state=42)
all_acc = []
all_auc = []

fold_preds = []
fold_labels = []
fold_preds_labels = []
for fold, (train_idx, val_idx) in enumerate(kf.split(all_index)):
    train_idx = [all_index[i] for i in train_idx]
    val_idx = [all_index[i] for i in val_idx]
    train_data = [all_data[i] for i in train_idx]
    train_label = [all_label[i] for i in train_idx]
    test_data = [all_data[i] for i in val_idx]
    test_label = [all_label[i] for i in val_idx]


    train_data = Data(train_data, train_label)
    test_data = Data(test_data, test_label)

    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=batch_size, shuffle=True)

    model = boundry_model(num_ep=3)
    loss_fun = nn.CrossEntropyLoss()
    opt = torch.optim.Adam(model.parameters(), lr=0.001)



    cur_fold_acc = -1
    cur_fold_auc = -1
    all_preds = []
    all_labels = []
    all_preds_labels = []
    for epoch in range(100):
        model.train()
        for batch in train_loader:
            pre_stage = batch[0]
            target_stage = batch[1]
            pre_level = batch[2]
            atact_value = batch[3]
            k4me3_value = batch[4]
            k9me3_value = batch[5]
            k27me3_value = batch[6]
            label = batch[7]

            ep_feature = torch.stack([k4me3_value, k9me3_value, k27me3_value], dim=0)
            ep_feature = ep_feature.T
            pre_stage = pre_stage
            target_stage = target_stage
            pre_level = pre_level

            pred_res = model(pre_stage, target_stage, pre_level, ep_feature.float())
            loss = loss_fun(pred_res, label)
            loss.backward()
            if math.isnan(loss.item()):
                print(ep_feature)
                exit(0)
            opt.step()
            opt.zero_grad()
        model.eval()
        cur_preds = []
        cur_labels = []
        cur_pred_labels = []
        for batch in test_loader:
            pre_stage = batch[0]
            target_stage = batch[1]
            pre_level = batch[2]
            atact_value = batch[3]
            k4me3_value = batch[4]
            k9me3_value = batch[5]
            k27me3_value = batch[6]
            label = batch[7]

            # ep_feature = torch.stack([k27me3_value], dim=0)
            ep_feature = torch.stack([k4me3_value, k9me3_value, k27me3_value], dim=0)
            # ep_feature = torch.stack([atact_value], dim=0)
            ep_feature = ep_feature.T
            pre_stage = pre_stage
            target_stage = target_stage
            pre_level = pre_level

            pred_res = model(pre_stage, target_stage, pre_level, ep_feature.float())
            pred_labels = torch.argmax(pred_res, dim=-1)
            pred_labels = pred_labels.numpy().tolist()
            cur_pred_labels.extend(pred_labels)
            bacth_acc = acc_metric(pred_res, label)
            batch_auc = auc_metric(pred_res, label)

            pred_res = pred_res.detach().numpy()[:, 1]
            pred_res = pred_res.tolist()
            cur_preds.extend(pred_res)
            label = label.detach().numpy().tolist()
            cur_labels.extend(label)

            fold_preds.extend(pred_res)
            fold_labels.extend(label)
            fold_preds_labels.extend(pred_labels)
        epoch_acc = acc_metric.compute()
        epoch_auc = auc_metric.compute()
        cur_fold_acc = max(epoch_acc, cur_fold_acc)
        cur_fold_auc = max(epoch_auc, cur_fold_auc)
        if cur_fold_auc > epoch_auc:
            all_preds = cur_preds
            all_labels = cur_labels
            all_preds_labels = cur_pred_labels
        acc_metric.reset()
        auc_metric.reset()
    all_acc.append(cur_fold_acc)
    all_auc.append(cur_fold_auc)
    fpr, tpr, thersholds = roc_curve(all_labels, all_preds, pos_label=1)
    roc_auc = auc(fpr, tpr)
    acc = accuracy_score(all_labels, all_preds_labels)
    f1 = f1_score(all_labels, all_preds_labels)
    plt.plot(fpr, tpr, 'k--', label='ROC (area = {0:.2f})'.format(roc_auc), lw=2)
    print('fold: {}, acc: {}, auc: {}, f1: {}'.format(fold, acc, roc_auc, f1))


    plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，以免和边缘重合，更好的观察图像的整体
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')  # 可以使用中文，但需要导入一些库即字体
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.show()
    plt.close()


fpr, tpr, thersholds = roc_curve(fold_labels, fold_preds, pos_label=1)
roc_auc = auc(fpr, tpr)
acc = accuracy_score(fold_labels, fold_preds_labels)
f1 = f1_score(fold_labels, fold_preds_labels)
plt.plot(fpr, tpr, 'k--', label='ROC (area = {0:.2f})'.format(roc_auc), lw=2)
plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，以免和边缘重合，更好的观察图像的整体
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')  # 可以使用中文，但需要导入一些库即字体
plt.title('ROC Curve')
plt.legend(loc="lower right")
plt.show()
plt.close()
print(roc_auc)
print(acc)
print(f1)
