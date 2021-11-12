import os
import numpy as np
import matplotlib.pyplot as plt
import random
import torch

def set_randomSeed(SEED=11):
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.backends.cudnn.deterministic = True

def check_folder(path):
    exist = os.path.exists(path)
    if not exist:
        os.makedirs(path)  # 若只用 mkdir()，如果父目錄不存在，會報錯
        print(path + "建立成功")
    else:
        print(path + "目錄已存在")
        
def plot(train_metirc, val_metric, metric_name, loss=False):
    plt.plot(train_metirc, label='train_%s' %metric_name)
    plt.plot(val_metric, label="val_%s" %metric_name)
    plt.xlabel('epochs')
    plt.ylabel(metric_name)
    if loss:
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='lower right')
    plt.show()