import torch
from torch import nn


class boundry_model(nn.Module):
    def __init__(self, num_ep):
        super(boundry_model, self).__init__()
        self.ep_feature_layer = nn.Linear(in_features=num_ep, out_features=16)
        self.stage_embed = nn.Embedding(embedding_dim=16, num_embeddings=8)
        self.level_embed = nn.Embedding(embedding_dim=16, num_embeddings=60)
        self.level_stage_project = nn.Linear(in_features=16, out_features=16)
        self.act = nn.Tanh()
        self.linear1 = nn.Linear(in_features=48, out_features=32)
        self.linear2 = nn.Linear(in_features=32, out_features=16)
        self.linear3 = nn.Linear(in_features=16, out_features=8)
        self.cls = nn.Linear(in_features=8, out_features=2)
        self.cls_act = nn.Softmax(dim=-1)
        # self.cls_act = nn.Sigmoid()


    def forward(self, pre_stage, target_stage, pre_level, ep_feature):
        target_stage_embedding = self.stage_embed(target_stage)
        pre_stage_level = pre_stage * 5 + pre_level
        pre_stage_level = pre_stage_level.int()
        level_stage_feature = self.level_embed(pre_stage_level)

        ep_feature = self.act(self.ep_feature_layer(ep_feature))

        level_stage_feature = self.act(self.level_stage_project(level_stage_feature))
        features = torch.cat((level_stage_feature, target_stage_embedding, ep_feature), dim=-1)
        # features =level_stage_feature
        features = self.act(self.linear1(features))
        features = self.act(self.linear2(features))
        features = self.act(self.linear3(features))
        cls = self.cls_act(self.cls(features))
        return cls


