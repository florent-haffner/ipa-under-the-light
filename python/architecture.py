from numpy.random import seed

import torch
from torch import Tensor

import torch.nn as nn
import torch.nn.functional as F

SEED_VALUE = 1
seed(SEED_VALUE)
torch.manual_seed(SEED_VALUE)

class BasicConv2d(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, **kwargs) -> None:
        super(BasicConv2d, self).__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, **kwargs)
        self.activation = nn.LeakyReLU(negative_slope=0.3)

    def forward(self, x: Tensor) -> Tensor:
        x = self.conv(x)
        x = self.activation(x)
        return x


class Module_35x35(nn.Module):
    def __init__(self, in_channels: int):
        super(Module_35x35, self).__init__()
        self.branch1 = nn.Sequential(
            nn.MaxPool1d(kernel_size=2),
            BasicConv2d(in_channels, in_channels*2, kernel_size=1, stride=2)
        )
        self.branch2 = nn.Sequential(
            BasicConv2d(in_channels, in_channels*2, kernel_size=1, stride=2),
            BasicConv2d(in_channels*2, in_channels*2, kernel_size=3, stride=1)
        )
        self.branch3 = nn.Sequential(
            BasicConv2d(in_channels, in_channels*2, kernel_size=1, stride=2),
            BasicConv2d(in_channels*2, in_channels*2, kernel_size=3, stride=2),
            BasicConv2d(in_channels*2, in_channels*2, kernel_size=3, stride=2)
        )
        self.branch4 = nn.Sequential(
            BasicConv2d(in_channels, in_channels*2, kernel_size=1, stride=2)
        )

    def forward(self, x: Tensor) -> Tensor:
        branch1 = self.branch1(x)
        branch2 = self.branch2(x)
        branch3 = self.branch3(x)
        branch4 = self.branch4(x)
        # print(f'1:{branch1.shape} - 2:{branch2.shape} - 3:{branch3.shape} - 4:{branch4.shape}')
        out = [branch1, branch2, branch3, branch4]
        out = torch.cat(out, 2)
        return out


class IPA(nn.Module):
    def __init__(self):
        super(IPA, self).__init__()
        self.stem = nn.Sequential(
            nn.Conv1d(1, 16, kernel_size=3, stride=2),
            nn.Conv1d(16, 16, kernel_size=3),
            nn.Conv1d(16, 32, kernel_size=3),
        )
        self.module_35x35 = Module_35x35(32)

        self.flatten = nn.Flatten()
        self.drop = nn.Dropout(.2)
        self.regressor = nn.LazyLinear(1) # The in_features argument of the Linear is inferred from the input.shape[-1].

    def forward(self, x):
        out = self.stem(x)
        out = self.module_35x35(out)
        out = self.flatten(out)
        out = self.drop(out)
        out = self.regressor(out)
        return out