import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# ============================================================
# 用户参数
# ============================================================

# 如果这个 py 文件放在 csv_out 文件夹里，就用 "."
# 如果 py 文件放在工程根目录，而 CSV 在 csv_out 里，就改成 "csv_out"
DATA_DIR = "."

# 可选：
# "phi"  看相场，最适合观察气泡
# "rho"  看密度场
# "T"    看温度场
FIELD = "phi"

NX = 160
NY = 100

# 动画播放速度
INTERVAL = 180          # 每帧间隔，单位 ms
FPS = 6                 # 保存 GIF 的帧率

# 是否保存 GIF
SAVE_GIF = True
GIF_NAME = "bubble_animation_phi.gif"

# 是否显示 phi=0 的气液界面线
DRAW_INTERFACE = True

# 你的最顶部 y=99 行此前可能是边界伪影，建议丢掉
DROP_TOP_ROW = True

# 只读取前多少帧；None 表示全部
MAX_FRAMES = None


# ============================================================
# 搜索 CSV 文件
# ============================================================

pattern = os.path.join(DATA_DIR, "fields_*.csv")
files = glob.glob(pattern)

if not files:
    raise FileNotFoundError(f"没有找到 CSV 文件：{pattern}")


def extract_step(filename):
    """
    从 fields_0000500.csv 中提取 500
    """
    name = os.path.basename(filename)
    match = re.search(r"fields_(\d+)\.csv", name)
    if match:
        return int(match.group(1))
    return -1


files = sorted(files, key=extract_step)

if MAX_FRAMES is not None:
    files = files[:MAX_FRAMES]

print(f"共找到 {len(files)} 个 CSV 场文件")
print("前几个文件：")
for f in files[:5]:
    print(" ", os.path.basename(f))


# ============================================================
# 读取单帧
# ============================================================

def load_frame(csv_file, field_name):
    """
    读取一个 fields_xxxxxxx.csv，
    返回：
        field_array: 需要显示的二维场
        phi_array: 相场二维数组，用于画 phi=0 界面
    """
    df = pd.read_csv(csv_file)

    required = {"x", "y", field_name}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_file} 缺少列：{missing}")

    # 去掉顶部伪边界行 y = NY - 1
    if DROP_TOP_ROW:
        df = df[df["y"] < NY - 1]

    ny_use = int(df["y"].max()) + 1

    field_array = np.full((ny_use, NX), np.nan, dtype=float)
    phi_array = None

    if "phi" in df.columns:
        phi_array = np.full((ny_use, NX), np.nan, dtype=float)

    for _, row in df.iterrows():
        x = int(row["x"])
        y = int(row["y"])

        if 0 <= x < NX and 0 <= y < ny_use:
            field_array[y, x] = row[field_name]

            if phi_array is not None:
                phi_array[y, x] = row["phi"]

    return field_array, phi_array


# ============================================================
# 读取第一帧
# ============================================================

first_field, first_phi = load_frame(files[0], FIELD)


# ============================================================
# 设置颜色范围
# ============================================================

if FIELD == "phi":
    vmin, vmax = -1.0, 1.0
    cmap = "coolwarm"

elif FIELD == "rho":
    all_min = np.inf
    all_max = -np.inf

    for f in files:
        arr, _ = load_frame(f, FIELD)
        all_min = min(all_min, np.nanmin(arr))
        all_max = max(all_max, np.nanmax(arr))

    vmin, vmax = all_min, all_max
    cmap = "viridis"

elif FIELD == "T":
    all_min = np.inf
    all_max = -np.inf

    for f in files:
        arr, _ = load_frame(f, FIELD)
        all_min = min(all_min, np.nanmin(arr))
        all_max = max(all_max, np.nanmax(arr))

    vmin, vmax = all_min, all_max
    cmap = "inferno"

else:
    all_min = np.inf
    all_max = -np.inf

    for f in files:
        arr, _ = load_frame(f, FIELD)
        all_min = min(all_min, np.nanmin(arr))
        all_max = max(all_max, np.nanmax(arr))

    vmin, vmax = all_min, all_max
    cmap = "viridis"


print(f"{FIELD} 色标范围：[{vmin:.6f}, {vmax:.6f}]")


# ============================================================
# 初始化画布
# ============================================================

fig, ax = plt.subplots(figsize=(9, 5.8))

image = ax.imshow(
    first_field,
    origin="lower",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
    interpolation="nearest"
)

colorbar = plt.colorbar(image, ax=ax)
colorbar.set_label(FIELD)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(0, NX - 1)
ax.set_ylim(0, NY - 2 if DROP_TOP_ROW else NY - 1)

title = ax.set_title("")

# 用字典保存 contour 对象，方便逐帧删除
contour_holder = {"obj": None}


# ============================================================
# 更新动画帧
# ============================================================

def update(frame_index):
    csv_file = files[frame_index]
    step = extract_step(csv_file)

    field_array, phi_array = load_frame(csv_file, FIELD)

    image.set_data(field_array)

    # 删除上一帧界面线
    if contour_holder["obj"] is not None:
        try:
            contour_holder["obj"].remove()
        except Exception:
            pass
        contour_holder["obj"] = None

    # 重新画当前帧 phi=0 等值线
    if DRAW_INTERFACE and phi_array is not None:
        try:
            contour_holder["obj"] = ax.contour(
                phi_array,
                levels=[0.0],
                colors="black",
                linewidths=1.2,
                origin="lower"
            )
        except Exception:
            contour_holder["obj"] = None

    title.set_text(f"{FIELD} field | step = {step}")

    return [image, title]


# ============================================================
# 生成动画
# ============================================================

animation = FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=INTERVAL,
    blit=False,
    repeat=True
)


# ============================================================
# 保存 GIF
# ============================================================

if SAVE_GIF:
    print(f"正在保存 GIF：{GIF_NAME}")
    writer = PillowWriter(fps=FPS)
    animation.save(GIF_NAME, writer=writer)
    print("GIF 保存完成")


plt.show()