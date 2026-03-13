import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
from pathlib import Path

# ======================== 可拖动标注类 ========================
class DraggableAnnotation:
    """可拖动的能垒标注"""
    def __init__(self, annotation, callback=None):
        self.annotation = annotation
        self.callback = callback
        self.draggable = annotation.draggable(True)
        # 监听拖动事件
        self.annotation.figure.canvas.mpl_connect('button_release_event', self._on_release)
        
    def _on_release(self, event):
        """鼠标释放时调用回调函数"""
        if self.callback is not None and event.inaxes == self.annotation.axes:
            self.callback()
    
    def on_key(self, event):
        """键盘方向键微调位置"""
        if event.key not in ['left', 'right', 'up', 'down']:
            return
        
        # 获取当前位置
        x, y = self.annotation.get_position()
        ax = self.annotation.axes
        
        # 计算步长 (坐标轴范围的1%)
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        dx = x_range * 0.01
        dy = y_range * 0.01
        
        # 根据按键调整位置
        if event.key == 'left':
            x -= dx
        elif event.key == 'right':
            x += dx
        elif event.key == 'up':
            y += dy
        elif event.key == 'down':
            y -= dy
        
        # 更新位置
        self.annotation.set_position((x, y))
        self.annotation.figure.canvas.draw_idle()
        
        if self.callback is not None:
            self.callback()

# ======================== 配置参数 ========================
# ---------- 数据设置 ----------
# 0K BSSW稳定结构
ENERGY_0K_ALLOY = -1818.4879      # 合金态能量 (eV)
ALPHA_0K_ALLOY = -0.171791         # 合金态alpha值

# 900K回火平均能量
ENERGY_900K_AVG = -1817.651325    # 偏析态能量 (eV)
ALPHA_900K_AVG = 0.125088          # 偏析态alpha值

# 从自由能曲线读取的数据 (alpha=-0.170 和 alpha=0.130)
# 300K数据 (已经是相对于合金态参考点的能垒)
ENERGY_300K_ALLOY_DEFAULT = 0.0      # 300K合金态参考点 (设为0)
ENERGY_300K_SEG_DEFAULT = 0.6005     # 300K能垒 (eV)
# 900K数据 (已经是相对于合金态参考点的能垒)
ENERGY_900K_ALLOY_DEFAULT = 0.0      # 900K合金态参考点 (设为0)
ENERGY_900K_SEG_DEFAULT = -0.4283    # 900K能垒 (eV)

# 这些值可以通过命令行参数覆盖
ENERGY_300K_ALLOY = None   
ENERGY_300K_SEG = None     
ENERGY_900K_ALLOY = None   
ENERGY_900K_SEG = None

# ---------- 绘图设置 ----------
# 颜色方案
COLOR_0K = '#1f77b4'       # 蓝色
COLOR_300K = '#ff7f0e'     # 橙色
COLOR_900K = '#d62728'     # 红色

# 台阶样式
STEP_LINE_WIDTH = 4        # 台阶线宽
STEP_LINE_STYLE = '-'      # 台阶线样式
HORIZONTAL_LINE_WIDTH = 4  # 水平线宽
VERTICAL_LINE_WIDTH = 2    # 垂直线宽
MARKER_SIZE = 100          # 端点标记大小
MARKER_STYLE = 'o'         # 标记样式

# 字体设置
FONT_FAMILY = 'Arial'
TICK_LABEL_SIZE = 28
AXIS_LABEL_SIZE = 34
LEGEND_SIZE = 26
TITLE_SIZE = 28

# 坐标轴设置
X_LABEL = 'State'
Y_LABEL = 'Relative Energy (eV)'
TICK_DIRECTION = 'out'
TICK_LENGTH = 6
TICK_WIDTH = 1.5
SPINE_WIDTH = 1.5

# 图片设置
FIGURE_WIDTH = 10
FIGURE_HEIGHT = 8
DPI = 300
TRANSPARENT_BG = True
OUTPUT_FORMAT = 'png'

# 图例设置
SHOW_LEGEND = True
LEGEND_FRAME = False       # 无图例框
LEGEND_LOCATION = 'upper right'

# 标题设置
SHOW_TITLE = False         # 不显示标题

# ==========================================================

# 解析命令行参数
parser = argparse.ArgumentParser(
    description='绘制自由能台阶图 (0K, 300K, 900K)',
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('--energy-300k-alloy', type=float, default=None,
                    help='300K合金态自由能 (eV, 从FES曲线获取)')
parser.add_argument('--energy-300k-seg', type=float, default=None,
                    help='300K偏析态自由能 (eV, 从FES曲线获取)')
parser.add_argument('--energy-900k-alloy', type=float, default=None,
                    help='900K合金态自由能 (eV, 从FES曲线获取)')
parser.add_argument('--energy-900k-seg', type=float, default=None,
                    help='900K偏析态自由能 (eV, 从FES曲线获取)')
parser.add_argument('--alpha-alloy', type=float, default=ALPHA_0K_ALLOY,
                    help=f'合金态alpha值 (默认: {ALPHA_0K_ALLOY})')
parser.add_argument('--alpha-seg', type=float, default=ALPHA_900K_AVG,
                    help=f'偏析态alpha值 (默认: {ALPHA_900K_AVG})')

args = parser.parse_args()

# 应用命令行参数
if args.energy_300k_alloy is not None:
    ENERGY_300K_ALLOY = args.energy_300k_alloy
else:
    ENERGY_300K_ALLOY = ENERGY_300K_ALLOY_DEFAULT
    
if args.energy_300k_seg is not None:
    ENERGY_300K_SEG = args.energy_300k_seg
else:
    ENERGY_300K_SEG = ENERGY_300K_SEG_DEFAULT
    
if args.energy_900k_alloy is not None:
    ENERGY_900K_ALLOY = args.energy_900k_alloy
else:
    ENERGY_900K_ALLOY = ENERGY_900K_ALLOY_DEFAULT
    
if args.energy_900k_seg is not None:
    ENERGY_900K_SEG = args.energy_900k_seg
else:
    ENERGY_900K_SEG = ENERGY_900K_SEG_DEFAULT
    
if args.alpha_alloy:
    ALPHA_0K_ALLOY = args.alpha_alloy
if args.alpha_seg:
    ALPHA_900K_AVG = args.alpha_seg

# 设置字体
mpl.rcParams['font.family'] = FONT_FAMILY
mpl.rcParams['font.sans-serif'] = [FONT_FAMILY]

# 计算相对能量 (以各自合金态为零点)
# 0K数据
E_0K_alloy_ref = ENERGY_0K_ALLOY  # 保存原始值用于打印
E_0K_alloy = 0.0  # 对齐为0点
E_0K_seg = ENERGY_900K_AVG - ENERGY_0K_ALLOY  # 相对于自己的合金态

# 300K数据
E_300K_alloy_ref = ENERGY_300K_ALLOY  # 保存原始值
E_300K_alloy = 0.0  # 对齐为0点
E_300K_seg = ENERGY_300K_SEG - ENERGY_300K_ALLOY  # 相对于自己的合金态

# 900K数据
E_900K_alloy_ref = ENERGY_900K_ALLOY  # 保存原始值
E_900K_alloy = 0.0  # 对齐为0点
E_900K_seg = ENERGY_900K_SEG - ENERGY_900K_ALLOY  # 相对于自己的合金态

# 准备数据
temperatures = ['0 K', '300 K', '900 K']
alpha_positions = [ALPHA_0K_ALLOY, ALPHA_900K_AVG]
alpha_labels = ['Alloy\n(α=-0.17)', 'Segregation\n(α=0.13)']

# 创建图形
fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

# 使用简单的X轴位置: 0=合金态, 1=偏析态
x_alloy = 0
x_seg = 1
step_width = 0.3  # 台阶的宽度

# 绘制0K台阶（黑色偏析态，灰色合金态）
# 合金态水平台阶（灰色）
ax.plot([x_alloy - step_width/2, x_alloy + step_width/2], [E_0K_alloy, E_0K_alloy], 
        color='#808080', linewidth=3, zorder=3)
# 偏析态水平台阶（黑色，带图例）
ax.plot([x_seg - step_width/2, x_seg + step_width/2], [E_0K_seg, E_0K_seg], 
        color='black', linewidth=3, zorder=3, label='0 K (DFT)')
# 连接虚线
ax.plot([x_alloy + step_width/2, x_seg - step_width/2], [E_0K_alloy, E_0K_seg], 
        color='black', linewidth=2, linestyle='--', alpha=0.5, zorder=2)
# 能垒标注（黑色，带eV单位，可拖动）
barrier_0K = E_0K_seg - E_0K_alloy
ann_0K = ax.annotate(f'{barrier_0K:.2f} eV', 
                     xy=(x_seg, E_0K_seg + 0.08),
                     ha='center', va='bottom', fontsize=24, color='black', zorder=5)
drag_0K = DraggableAnnotation(ann_0K)

# 绘制300K台阶（蓝色偏析态,灰色合金态）
# 合金态水平台阶（灰色）
ax.plot([x_alloy - step_width/2, x_alloy + step_width/2], [E_300K_alloy, E_300K_alloy], 
        color='#808080', linewidth=3, zorder=3)
# 偏析态水平台阶（蓝色，带图例）
ax.plot([x_seg - step_width/2, x_seg + step_width/2], [E_300K_seg, E_300K_seg], 
        color='#1f77b4', linewidth=3, zorder=3, label='300 K (MD)')
# 连接虚线
ax.plot([x_alloy + step_width/2, x_seg - step_width/2], [E_300K_alloy, E_300K_seg], 
        color='#1f77b4', linewidth=2, linestyle='--', alpha=0.5, zorder=2)
# 能垒标注（蓝色，带eV单位，可拖动）
barrier_300K = E_300K_seg - E_300K_alloy
ann_300K = ax.annotate(f'{barrier_300K:.2f} eV', 
                       xy=(x_seg, E_300K_seg + 0.08),
                       ha='center', va='bottom', fontsize=24, color='#1f77b4', zorder=5)
drag_300K = DraggableAnnotation(ann_300K)

# 绘制900K台阶（红色偏析态,灰色合金态）
# 合金态水平台阶（灰色）
ax.plot([x_alloy - step_width/2, x_alloy + step_width/2], [E_900K_alloy, E_900K_alloy], 
        color='#808080', linewidth=3, zorder=3)
# 偏析态水平台阶（红色，带图例）
ax.plot([x_seg - step_width/2, x_seg + step_width/2], [E_900K_seg, E_900K_seg], 
        color='#d62728', linewidth=3, zorder=3, label='900 K (MD)')
# 连接虚线
ax.plot([x_alloy + step_width/2, x_seg - step_width/2], [E_900K_alloy, E_900K_seg], 
        color='#d62728', linewidth=2, linestyle='--', alpha=0.5, zorder=2)
# 能垒标注（红色，带eV单位，可拖动）
barrier_900K = E_900K_seg - E_900K_alloy
# 如果是负值，标注在上方；否则标注在下方
if barrier_900K < 0:
    ann_900K = ax.annotate(f'{barrier_900K:.2f} eV', 
                           xy=(x_seg, E_900K_seg + 0.08),
                           ha='center', va='bottom', fontsize=24, color='#d62728', zorder=5)
else:
    ann_900K = ax.annotate(f'{barrier_900K:.2f} eV', 
                           xy=(x_seg, E_900K_seg + 0.08),
                           ha='center', va='bottom', fontsize=24, color='#d62728', zorder=5)
drag_900K = DraggableAnnotation(ann_900K)

# 设置坐标轴标签（不加粗）
ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel(Y_LABEL, fontsize=AXIS_LABEL_SIZE)

# 设置刻度
ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE,
               direction=TICK_DIRECTION, length=TICK_LENGTH, width=TICK_WIDTH)

# 设置X轴范围和刻度
x_margin = 0.3
ax.set_xlim(x_alloy - x_margin, x_seg + x_margin)
ax.set_xticks([x_alloy, x_seg])
ax.set_xticklabels(['Alloy', 'Segregation'])

# 设置Y轴范围和刻度（对称整数刻度）
y_values = [E_0K_alloy, E_0K_seg, E_300K_alloy, E_300K_seg, E_900K_alloy, E_900K_seg]
y_min = min(y_values)
y_max = max(y_values)

# 使用对称的整数刻度
ax.set_yticks([-0.5, 0.0, 0.5, 1.0])
ax.set_ylim(-0.6, 1.1)

# 设置边框（4个框都显示）
for spine in ax.spines.values():
    spine.set_linewidth(SPINE_WIDTH)
    spine.set_visible(True)

# 关闭网格线
ax.grid(False)

# 图例
if SHOW_LEGEND:
    legend = ax.legend(fontsize=LEGEND_SIZE, frameon=LEGEND_FRAME, 
                      loc=LEGEND_LOCATION, framealpha=0.9)
    legend.set_draggable(True)

# 添加键盘事件监听（支持方向键微调）
def on_key(event):
    """键盘事件处理"""
    drag_0K.on_key(event)
    drag_300K.on_key(event)
    drag_900K.on_key(event)

fig.canvas.mpl_connect('key_press_event', on_key)

# 不添加网格线
# ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, zorder=1)

# 设置背景透明
if TRANSPARENT_BG:
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

# 紧凑布局
plt.tight_layout()

# 准备保存路径
output_dir = Path(r'C:\Users\11207\OneDrive\02_Code\work1-PtSnO\workflow\results\part5')
output_dir.mkdir(parents=True, exist_ok=True)
output_filename = f'energy_barrier_steps.{OUTPUT_FORMAT}'
output_path = output_dir / output_filename

save_params = {
    'dpi': DPI,
    'bbox_inches': 'tight',
}
if TRANSPARENT_BG:
    save_params['facecolor'] = 'none'
    save_params['edgecolor'] = 'none'
    save_params['transparent'] = True

# 定义窗口关闭时的保存函数
def on_close(event):
    """窗口关闭时保存图片"""
    fig.savefig(output_path, **save_params)
    print("\n" + "=" * 70)
    print("图片已保存!")
    print("=" * 70)
    print(f"输出路径: {output_path}")
    print("-" * 70)
    print("能量数据 (合金态对齐为0点):")
    print(f"  0K 合金态: {E_0K_alloy:.2f} eV (原始: {E_0K_alloy_ref:.2f} eV)")
    print(f"  0K 偏析态: {E_0K_seg:.2f} eV")
    print(f"  0K能垒: ΔE = {E_0K_seg - E_0K_alloy:.2f} eV")
    print("-" * 70)
    print(f"  300K 合金态: {E_300K_alloy:.2f} eV (原始: {E_300K_alloy_ref:.2f} eV)")
    print(f"  300K 偏析态: {E_300K_seg:.2f} eV")
    print(f"  300K能垒: ΔF = {E_300K_seg - E_300K_alloy:.2f} eV")
    print("-" * 70)
    print(f"  900K 合金态: {E_900K_alloy:.2f} eV (原始: {E_900K_alloy_ref:.2f} eV)")
    print(f"  900K 偏析态: {E_900K_seg:.2f} eV")
    print(f"  900K能垒: ΔF = {E_900K_seg - E_900K_alloy:.2f} eV")
    print("=" * 70)

# 绑定窗口关闭事件
fig.canvas.mpl_connect('close_event', on_close)

# 显示图形
print("=" * 70)
print("能量台阶图窗口已打开!")
print("=" * 70)
print("提示: 可以拖动图例到合适位置")
print("      可以拖动能垒数值标注到合适位置")
print("      支持键盘方向键微调标注位置 (↑↓←→)")
print("      关闭窗口后将自动保存图片")
print("-" * 70)
print("当前数据 (所有合金态对齐为0点):")
print(f"  0K DFT: 能垒 ΔE = {E_0K_seg:.2f} eV")
print(f"  300K MD: 能垒 ΔF = {E_300K_seg:.2f} eV")
print(f"  900K MD: 能垒 ΔF = {E_900K_seg:.2f} eV")
print("=" * 70)

plt.show()
