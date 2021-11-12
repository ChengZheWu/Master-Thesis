import torch.nn as nn
import torch

class RSKblock(nn.Module):
    def __init__(self, inplanes, planes, stride=1, M=2, r=8):
        super(RSKblock,self).__init__()
        self.M = M
        self.planes = planes
        midplanes = planes // 2
        self.convs = nn.ModuleList()
        for i in range(M):
            self.convs.append(nn.Sequential(
                nn.Conv2d(inplanes, midplanes, kernel_size=1, stride=1, bias=False),
                nn.BatchNorm2d(midplanes),
                nn.ReLU(inplace=True),
                nn.Conv2d(midplanes, midplanes, kernel_size=3+2*i, stride=stride, padding=1+i, bias=False),
                nn.BatchNorm2d(midplanes),
                nn.ReLU(inplace=True),
                nn.Conv2d(midplanes, planes, kernel_size=1, stride=1, bias=False),
                nn.BatchNorm2d(planes),
            ))
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.fc1 =nn.Sequential(
            nn.Linear(planes, planes//r),
        )
        self.fc2 = nn.ModuleList([])
        for i in range(M+1):
            self.fc2.append(nn.Linear(planes//r, planes))
        self.softmax = nn.Softmax(dim=1)
        self.reluf = nn.ReLU(inplace=False)
        self.relu = nn.ReLU(inplace=True)
        
        if stride == 1 and inplanes == planes:
            self.downsample = None
        else:
            self.downsample = nn.Sequential(
                nn.Conv2d(inplanes, planes, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(planes),
            )
        
    def forward(self, x):
        batch_size = x.size(0)
        
        # Split
        identity = x
        if self.downsample is not None:
            identity = self.downsample(x)
        output = identity.unsqueeze(dim=1)
        for conv in self.convs:
            output = torch.cat((output, conv(x).unsqueeze(dim=1)), dim=1)
        r_output = self.reluf(output)
        
        # Fuse
        U = torch.sum(r_output, dim=1)  # element-wise summation
    
        s = self.avgpool(U)
        s = s.view(batch_size, -1)
        z = self.fc1(s)
        for i, fc in enumerate(self.fc2):
            if i == 0:
                abc = fc(z).unsqueeze_(dim=1)
            else:
                abc = torch.cat([abc, fc(z).unsqueeze_(dim=1)], dim=1)
        abc = self.softmax(abc)
        abc = abc.view(batch_size, self.M+1, self.planes, 1, 1)
        
        # Select
        V = torch.sum(output*abc, dim=1)
        return self.relu(V)

class MRSKNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.block = nn.Sequential(
            RSKblock(1, 32, 2),
            RSKblock(32, 32),
            nn.Dropout2d(0.1),
            RSKblock(32, 64, 2),
            RSKblock(64, 64),
            nn.Dropout2d(0.1),
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