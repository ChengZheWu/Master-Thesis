import torch.nn as nn
import torch


class BasicBlock(nn.Module):
    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(BasicBlock, self).__init__()
        self.conv1 = nn.Conv2d(inplanes, planes, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(planes, planes, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(planes)
        if stride == 1 and inplanes == planes:
            self.downsample = None
        else:
            self.downsample = nn.Sequential(
                nn.Conv2d(inplanes, planes, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(planes),
            )

    def forward(self, x):
        identity = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            identity = self.downsample(x)

        out += identity
        out = self.relu(out)

        return out
    
class Bottleneck(nn.Module):
    def __init__(self, inplanes, planes, stride=1):
        super(Bottleneck, self).__init__()
        midplanes = planes // 2
        self.conv = nn.Sequential(
            nn.Conv2d(inplanes, midplanes, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm2d(midplanes),
            nn.ReLU(inplace=True),
            nn.Conv2d(midplanes, midplanes, kernel_size=3, stride=stride, padding=1, bias=False),
            nn.BatchNorm2d(midplanes),
            nn.ReLU(inplace=True),
            nn.Conv2d(midplanes, planes, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm2d(planes)
        )
        self.relu = nn.ReLU(inplace=True)
        if stride == 1 and inplanes == planes:
            self.downsample = None
        else:
            self.downsample = nn.Sequential(
                nn.Conv2d(inplanes, planes, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(planes),
            )
    
    def forward(self, x):
        
        identity = x
        if self.downsample is not None:
            identity = self.downsample(x)
        out = self.conv(x)
        out = out + identity
        return self.relu(out)


# class BasicResNet(nn.Module):
#     def __init__(self):
#         super().__init__()
#         self.block = nn.Sequential(
#             BasicBlock(2, 32, 2),
#             BasicBlock(32, 32),
#             nn.Dropout2d(0.2),
#             BasicBlock(32, 64, 2),
#             BasicBlock(64, 64),
#             nn.Dropout2d(0.2),
#         )
#         self.maxpool = nn.AdaptiveMaxPool2d(1)
#         self.fc = nn.Linear(192, 1)
    
#     def forward(self, x):
        
#         out1 = self.block(x[:, 0])
#         out1 = self.maxpool(out1)
#         out2 = self.block(x[:, 1])
#         out2 = self.maxpool(out2)
#         out3 = self.block(x[:, 2])
#         out3 = self.maxpool(out3)
        
#         out = torch.cat((out1, out2, out3), 1)
#         out = out.view(out.size(0), -1)
#         out = self.fc(out)
#         out = torch.sigmoid(out)
#         return out
    
class BasicResNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.block = nn.Sequential(
            Bottleneck(2, 32, 2),
            Bottleneck(32, 32),
            nn.Dropout2d(0.2),
            Bottleneck(32, 64, 2),
            Bottleneck(64, 64),
            nn.Dropout2d(0.2),
        )
        self.maxpool = nn.AdaptiveMaxPool2d(1)
        self.fc = nn.Linear(192, 1)
    
    def forward(self, x):
        
        out1 = self.block(x[:, 0])
        out1 = self.maxpool(out1)
        out2 = self.block(x[:, 1])
        out2 = self.maxpool(out2)
        out3 = self.block(x[:, 2])
        out3 = self.maxpool(out3)
        
        out = torch.cat((out1, out2, out3), 1)
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        out = torch.sigmoid(out)
        return out