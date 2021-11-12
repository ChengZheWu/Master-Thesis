import torch.nn as nn
import torch


def conv3x3(in_planes, out_planes, stride=1):
    return nn.Conv2d(in_planes, out_planes, kernel_size=3, stride=stride, padding=1, bias=False)

class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(BasicBlock, self).__init__()
        self.conv1 = conv3x3(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(planes, planes)
        self.bn2 = nn.BatchNorm2d(planes)
        self.downsample = downsample
        self.stride = stride

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
    

# Copied from https://github.com/heykeetae/Self-Attention-GAN
class SelfAttn(nn.Module):
    """ Self attention Layer"""

    def __init__(self, in_dim):
        super().__init__()
        self.chanel_in = in_dim

        self.query_conv = nn.Conv2d(in_channels=in_dim, out_channels=in_dim // 8, kernel_size=1)
        self.key_conv = nn.Conv2d(in_channels=in_dim, out_channels=in_dim // 8, kernel_size=1)
        self.value_conv = nn.Conv2d(in_channels=in_dim, out_channels=in_dim, kernel_size=1)
        self.gamma = nn.Parameter(torch.zeros(1))

        self.softmax = nn.Softmax(dim=-1)  #

    def forward(self, x):
        """
            inputs :
                x : input feature maps( B X C X W X H)
            returns :
                out : self attention value + input feature
                attention: B X N X N (N is Width*Height)
        """
        m_batchsize, C, width, height = x.size()
        proj_query = self.query_conv(x).view(m_batchsize, -1, width * height).permute(0, 2, 1)  # B X CX(N)
        proj_key = self.key_conv(x).view(m_batchsize, -1, width * height)  # B X C x (*W*H)
        energy = torch.bmm(proj_query, proj_key)  # transpose check
        attention = self.softmax(energy)  # BX (N) X (N)
        proj_value = self.value_conv(x).view(m_batchsize, -1, width * height)  # B X C X N

        out = torch.bmm(proj_value, attention.permute(0, 2, 1))
        out = out.view(m_batchsize, C, width, height)

        out = self.gamma * out + x


        return out

class LocalGlobalNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.resnetAttn = nn.Sequential(
            BasicBlock(1, 32),
            SelfAttn(32),
            nn.Dropout(0.1),
            BasicBlock(32, 32),
            SelfAttn(32),
            nn.Dropout(0.2),
        )
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.fc1 = nn.Linear(96, 1)


    def forward(self, x):

        x1 = self.resnetAttn(x[:, 0])     
        x1 = self.avgpool(x1)
        x2 = self.resnetAttn(x[:, 1])
        x2 = self.avgpool(x2)
        x3 = self.resnetAttn(x[:, 2])
        x3 = self.avgpool(x3)
        
        out = torch.cat((x1, x2, x3), 1)
#         out = self.avgpool(out)
        out = out.view(out.size(0), -1)
        out = self.fc1(out)
        out = torch.sigmoid(out)
        return out