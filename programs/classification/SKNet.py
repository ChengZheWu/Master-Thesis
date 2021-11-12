import torch.nn as nn
import torch

class SKconv(nn.Module):
    def __init__(self, inplanes, planes, stride=1, M=2, r=8):
        super(SKconv,self).__init__()
        self.M = M
        self.planes = planes
        self.convs = nn.ModuleList()
        for i in range(M):
            self.convs.append(nn.Sequential(
                nn.Conv2d(inplanes, planes, kernel_size=3+2*i, stride=stride, padding=1+i, bias=False),
                nn.BatchNorm2d(planes),
                nn.ReLU(inplace=True),
            ))
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.fc1 =nn.Sequential(
                nn.Linear(planes, planes//r))
        self.fc2 = nn.ModuleList([])
        for i in range(M):
            self.fc2.append(
                nn.Linear(planes//r, planes)
            )
        self.softmax = nn.Softmax(dim=1)
        self.relu = nn.ReLU(inplace=True)
        
    def forward(self, x):
        batch_size = x.size(0)
        
        # Split
        for i, conv in enumerate(self.convs):
            if i == 0:
                output = conv(x).unsqueeze(dim=1)
            else:
                output = torch.cat((output, conv(x).unsqueeze(dim=1)), dim=1)
        
        # Fuse
        U = torch.sum(output, dim=1)  # element-wise summation
    
        s = self.avgpool(U)
        s = s.view(batch_size, -1)
        z = self.fc1(s)
        for i, fc in enumerate(self.fc2):
            if i == 0:
                ab = fc(z).unsqueeze_(dim=1)
            else:
                ab = torch.cat([ab, fc(z).unsqueeze_(dim=1)], dim=1)
        ab = self.softmax(ab)
        ab = ab.view(batch_size, self.M, self.planes, 1, 1)

        # Select
        V = torch.sum(output*ab, dim=1)
        return V
    
class SKblock(nn.Module):
    def __init__(self, inplanes, planes, stride=1):
        super(SKblock, self).__init__()
        midplanes = planes // 2
        self.conv = nn.Sequential(
            nn.Conv2d(inplanes, midplanes, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm2d(midplanes),
            nn.ReLU(inplace=True),
            SKconv(midplanes, midplanes, stride),
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

class BasicSKNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.block = nn.Sequential(
            SKblock(2, 32, 2),
            SKblock(32, 32),
            nn.Dropout2d(0.2),
            SKblock(32, 64, 2),
            SKblock(64, 64),
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