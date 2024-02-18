# from pandas import DataFrame
# from scipy.stats import uniform
# from scipy.stats import randint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import colorsys
import random

def get_n_hls_colors(num):
    hls_colors = []
    i = 0
    step = 360.0 / num
    while i < 360:
        h = i
        s = 90 + random.random() * 10
        l = 50 + random.random() * 10
        _hlsc = [h / 360.0, l / 100.0, s / 100.0]
        hls_colors.append(_hlsc)
        i += step

    return hls_colors

def ncolors(num):
    rgb_colors = []
    if num < 1:
        return rgb_colors
    hls_colors = get_n_hls_colors(num)
    for hlsc in hls_colors:
        _r, _g, _b = colorsys.hls_to_rgb(hlsc[0], hlsc[1], hlsc[2])
        r, g, b = [int(x * 255.0) for x in (_r, _g, _b)]
        rgb_colors.append([r, g, b])

    return rgb_colors

def color(value):
    digit = list(map(str, range(10))) + list("ABCDEF")
    if isinstance(value, tuple):
        string = '#'
        for i in value:
            a1 = i // 16
            a2 = i % 16
            string += digit[a1] + digit[a2]
        return string
    elif isinstance(value, str):
        a1 = digit.index(value[1]) * 16 + digit.index(value[2])
        a2 = digit.index(value[3]) * 16 + digit.index(value[4])
        a3 = digit.index(value[5]) * 16 + digit.index(value[6])
        return (a1, a2, a3)

color_map = list(map(lambda x: color(tuple(x)), ncolors(6)))
print(color_map)
weight_index = 0
for column in ['A', 'B', 'C', 'D', 'E', 'F']:
    filename = r"result_r\SVM_weight_lc0303.xls"
    weight = pd.read_excel(filename, header=None, usecols=column)
    y = weight.values.tolist()

    x = []
    for num in range(400):
        x.append(math.floor(num / 10) + 1 - 0.4 + 0.05*weight_index)
    x = np.array(x)

    weight_tran = []
    for ywei in y:
        weight_tran.append(ywei[0])
    weight_tran = np.array(weight_tran)
    weight_tran = abs(weight_tran)

    color_value = []
    for color_num in range(400):
        color_value.append(color_map[weight_index])
    color_value = np.array(color_value)

    plt.scatter(x, weight_tran, c=color_value, alpha=1, cmap=plt.cm.rainbow, edgecolors='none', s=6, label='{}'.format(weight_index+1))
    # plt.legend(loc=(1.02, 0.58), ncol=1, fontsize=12, markerscale=2, handlelength=1.5, handletextpad=0.2, columnspacing=1.5)# frameon=False,
    weight_index = weight_index + 1

    # x = []
    # for num in range(7200):
    #     x.append(math.floor((num % 400)/10) + 1)
    # x = np.array(x)
    #
    # weight_tran = []
    # for ywei_num in range(18):
    #     for ywei in y:
    #         weight_tran.append(ywei[ywei_num])
    # weight_tran = np.array(weight_tran)
    # weight_tran = abs(weight_tran)
    #
    # color_value = []
    # for color_num in range(7200):
    #     color_value.append(math.floor(color_num/400) + 1)
    # color_value = np.array(color_value)
    #
    # plt.scatter(x, weight_tran, c=color_value, cmap=plt.cm.rainbow, edgecolors='none', s=3)
    #
    # weight_index = weight_index+1

plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 显示汉字
# plt.xlabel('node')  # x轴标题
# plt.ylabel('weight')  # y轴标题
# plt.title('the SVM weight in different nodes scale1')  # 折线图标题
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')
plt.ylim((-0.03, 1.0))
# plt.legend(['multi-scale', 'single-scale', 'global', 'original'], loc=(1.1, 4))  # 设置折线名称
plt.subplots_adjust(hspace=0.3, wspace=0.3, right=0.8)
plt.savefig('weight_node_scale1_220303.svg', bbox_inches='tight')
plt.show()  # 显示折线图
