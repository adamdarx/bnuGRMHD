import numpy as np
import sys

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MaxNLocator

N1, N2, N3, NPRIM = 128, 128, 4, 10
a = 0.9375
R0 = 0
h = 0

rhor = 1 + np.sqrt(1-a**2)

small = 1e-16
PI = 3.14159265358979323846

X1min = 0.19325057145871735
X1max = 7.824046010856292
X2min = 0
X2max = 2 * PI
X3min = 0
X3max = 2 * PI

dx1 = (X1max - X1min) / N1
dx2 = (X2max - X2min) / N2
dx3 = (X3max - X3min) / N3

#MKS grid
X1, X2, X3 = np.meshgrid(np.linspace(X1min, X1max, N1), np.linspace(X2min, X2max, N2), np.linspace(X3min, X3max, N3), indexing='ij')
#KS grid
r = np.exp(X1) + R0
theta = X2 + h / 2 * np.sin(2 * X2)
phi = X3

i_hor = np.where(r[:, 0, 0] >= rhor)[0][0]

def read_bin(file_path):
    total_elements = N1 * N2 * N3 * NPRIM
    data = np.fromfile(file_path, dtype=np.float64, count=total_elements)
    data = data.reshape((N1, N2, N3, NPRIM))
    return data

#black hole
def black_hole():
    x_b = np.linspace(-rhor,rhor,100)
    z1_b = np.sqrt(rhor**2 - x_b**2)
    z2_b = -z1_b
    return x_b, z1_b, z2_b

def rho_xz(lim):
    #get data lim from plot lim
    Rout = 1.5 * lim
    i_out = np.where(r[:, 0, 0] >= Rout)[0][0]
    #get plot grid
    x_p = r[i_hor:i_out,:,0] * np.sin(theta[i_hor:i_out,:,0])
    z_p = r[i_hor:i_out,:,0] * np.cos(theta[i_hor:i_out,:,0])
    #rho: average of phi
    rho_p = rho[i_hor:i_out].sum(2) / N3
    #black hole
    x_b, z1_b, z2_b = black_hole()
    #plot
    plt.clf()
    fig, ax = plt.subplots(figsize=(20,20), dpi=100)
    #contour of rho
    level = np.linspace(-1, 0, 100)
    contour = ax.contourf(x_p, z_p, np.log10(rho_p), levels=level, cmap='viridis', extend='both')
    #horizon
    ax.fill_between(x_b, z1_b, z2_b, color='black')
    #plot range
    ax.set_xlim(0, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal','box')
    #title and label
    size = 30
    ax.set_title(r'$t$=%03d[$r_{\mathrm{g}}/c$]' %t, fontsize=size)
    ax.set_xlabel(r'$x[r_{\mathrm{g}}]$', fontsize=size)
    ax.set_ylabel(r'$z[r_{\mathrm{g}}]$    ', rotation=90, fontsize=size)
    ax.tick_params(labelsize = size)
    #colorbar
    cb = fig.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.set_ylabel(r'log$_{10} \rho$', fontsize=size, rotation=90, labelpad=20)
    cb.formatter = FormatStrFormatter('%.1f')
    cb.ax.tick_params(labelsize = size)
    cb.update_ticks()
    plt.savefig('./pic/rho_xz%04d.png' %t, bbox_inches='tight')
    print(t)
    plt.cla()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        for t in range(int(sys.argv[1]) + 1):
            file_path = f"./data/data{t:0>4d}.bin"
            data = read_bin(file_path)
            has_nan = np.any(np.isnan(data))
            if has_nan:
                print("数组中存在NaN值")
                nan_indices = np.where(np.isnan(data))
                print("NaN 的索引位置 (三维数组)：")
                for i in range(len(nan_indices[0])):
                    print(f"NaN 位于索引 ({nan_indices[0][i]}, {nan_indices[1][i]}, {nan_indices[2][i]}, {nan_indices[3][i]})")
            else:
                print("数组中不存在NaN值")
            rho = data[:, :, :, 0]
            ug = data[:, :, :, 1]
            uu = data[:, :, :, 2:5+1]
            B = data[:, :, :, 6:8+1]
            bsq = data[:, :, :, 9]
            rho_xz(20)
    else:
        for t in range(int(sys.argv[1]) + 1, int(sys.argv[2]) + 1):
            file_path = f"./data/data{t:0>4d}.bin"
            data = read_bin(file_path)
            has_nan = np.any(np.isnan(data))
            if has_nan:
                print("数组中存在NaN值")
                nan_indices = np.where(np.isnan(data))
                print("NaN 的索引位置 (三维数组)：")
                for i in range(len(nan_indices[0])):
                    print(f"NaN 位于索引 ({nan_indices[0][i]}, {nan_indices[1][i]}, {nan_indices[2][i]}, {nan_indices[3][i]})")
            else:
                print("数组中不存在NaN值")
            rho = data[:, :, :, 0]
            ug = data[:, :, :, 1]
            uu = data[:, :, :, 2:5+1]
            B = data[:, :, :, 6:8+1]
            bsq = data[:, :, :, 9]
            rho_xz(20)




