import torch.nn as nn
import torch.nn.functional as F


def conv(in_planes, out_planes, stride=1):
    """3x3 convolution with padding"""
    return nn.Conv2d(in_planes, out_planes, kernel_size=3, padding=1, stride=stride, bias=False)


class BasicConv2d(nn.Module):

    def __init__(self, in_channels, out_channels, **kwargs):
        super(BasicConv2d, self).__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, bias=False, **kwargs)
        self.bn = nn.BatchNorm2d(out_channels, eps=0.001)

    def forward(self, x):
        x = self.conv(x)
        x = self.bn(x)
        return F.relu(x, inplace=True)


class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(BasicBlock, self).__init__()
        self.conv1 = conv(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv(planes, planes)
        self.bn2 = nn.BatchNorm2d(planes)
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.conv2(out)
        out = self.bn2(out)
        out += residual
        out = self.relu(out)

        return out


class ResNet(nn.Module):

    def __init__(self, in_channels, block, layers):
        self.inplanes = 1
        super(ResNet, self).__init__()
        self.Context_Conv2d_0a = BasicConv2d(in_channels, 2 * in_channels, kernel_size=(1, 3), stride=(1, 1))
        self.Context_Conv2d_1a = BasicConv2d(2 * in_channels, 1, kernel_size=(1, 3), stride=(1, 2))

    def forward(self, x):
        x = self.Context_Conv2d_0a(x)
        x = self.Context_Conv2d_1a(x)

        return x


def resnet18_custom(input_channels):
    """Constructs a ResNet-18 model.
    """
    model = ResNet(input_channels, BasicBlock, [2])

    return model
