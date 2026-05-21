import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 用户参数
# ============================================================

DATA_DIR = "."      # 如果脚本在 csv_out 里，改成 "."
OUTPUT_DIR = "bubble_metrics_out"

NX = 160
NY = 100

# 判断气泡区域的相场阈值
# phi < 0 认为是气相
PHI_THRESHOLD = 0.0

# 壁面附近统计高度
# y <= WALL_Y_MAX 的气相格子认为是壁面气相/未完全再润湿区域
WALL_Y_MAX = 3

# 格子尺寸。如果你没有物理单位，就保持 1.0
DX = 1.0
DY = 1.0

# 是否去掉顶部固壁/边界行
DROP_TOP_ROW = True


# ============================================================
# 查找 CSV 文件
# ============================================================

def extract_step(filename):
    name = os.path.basename(filename)
    m = re.search(r"fields_(\d+)\.csv", name)
    if m:
        return int(m.group(1))
    return -1


def find_csv_files():
    pattern = os.path.join(DATA_DIR, "fields_*.csv")
    files = sorted(glob.glob(pattern), key=extract_step)

    if not files:
        raise FileNotFoundError(f"没有找到 CSV 文件：{pattern}")

    print(f"找到 {len(files)} 个 CSV 文件")
    print("前几个文件：")
    for f in files[:5]:
        print(" ", os.path.basename(f))

    return files


# ============================================================
# 计算单帧气泡指标
# ============================================================

def compute_metrics_from_file(csv_file):
    df = pd.read_csv(csv_file)

    required_cols = {"x", "y", "phi"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{csv_file} 缺少列：{missing}")

    if DROP_TOP_ROW:
        df = df[df["y"] < NY - 1]

    step = extract_step(csv_file)

    vapor = df[df["phi"] < PHI_THRESHOLD].copy()

    vapor_cells = len(vapor)

    if vapor_cells > 0:
        # 面积，二维中就是气相格子数乘以 dx*dy
        area = vapor_cells * DX * DY

        # 质心
        xc = float((vapor["x"] * DX).mean())
        yc = float((vapor["y"] * DY).mean())

        # 2D 等效直径
        # A = pi D^2 / 4 -> D = sqrt(4A/pi)
        d_eq = np.sqrt(4.0 * area / np.pi)

        # 最高和最低气泡位置
        y_min = float(vapor["y"].min())
        y_max = float(vapor["y"].max())
        x_min = float(vapor["x"].min())
        x_max = float(vapor["x"].max())
    else:
        area = 0.0
        xc = np.nan
        yc = np.nan
        d_eq = 0.0
        y_min = np.nan
        y_max = np.nan
        x_min = np.nan
        x_max = np.nan

    wall_vapor_cells = len(df[(df["phi"] < PHI_THRESHOLD) & (df["y"] <= WALL_Y_MAX)])

    return {
        "step": step,
        "vapor_cells": vapor_cells,
        "area": area,
        "x_centroid": xc,
        "y_centroid": yc,
        "D_eq": d_eq,
        "wall_vapor_cells": wall_vapor_cells,
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
    }


# ============================================================
# 简单识别脱落时刻
# ============================================================

def estimate_departure_step(metrics):
    """
    一个简单判据：
    wall_vapor_cells 明显下降，并且 y_centroid 开始上升。
    这个不是严格物理定义，只是辅助判断。
    """

    df = metrics.copy()

    if len(df) < 5:
        return None

    wall = df["wall_vapor_cells"].to_numpy()
    yc = df["y_centroid"].to_numpy()
    steps = df["step"].to_numpy()

    # 忽略没有气泡的帧
    valid = np.isfinite(yc) & (df["vapor_cells"].to_numpy() > 0)
    if valid.sum() < 5:
        return None

    # 壁面气相峰值
    wall_max = np.nanmax(wall)
    if wall_max <= 0:
        return None

    # 条件：壁面气相降到峰值的 40% 以下，同时质心高于早期质心
    early_yc = np.nanmedian(yc[valid][:max(3, valid.sum() // 5)])

    for i in range(1, len(df)):
        if not np.isfinite(yc[i]):
            continue

        wall_drop = wall[i] < 0.40 * wall_max
        centroid_lift = yc[i] > early_yc + 5.0

        if wall_drop and centroid_lift:
            return int(steps[i])

    return None


# ============================================================
# 画图函数
# ============================================================

def save_single_plot(df, x_col, y_col, ylabel, title, filename, departure_step=None):
    plt.figure(figsize=(8, 5))

    plt.plot(df[x_col], df[y_col], linewidth=2)

    if departure_step is not None:
        plt.axvline(departure_step, linestyle="--", linewidth=1.5)
        plt.text(
            departure_step,
            np.nanmax(df[y_col]) * 0.95 if np.nanmax(df[y_col]) > 0 else 0.0,
            "departure",
            rotation=90,
            va="top"
        )

    plt.xlabel("step")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(path, dpi=200)
    plt.close()

    print(f"保存：{path}")


def save_4in1_plot(df, departure_step=None):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    plots = [
        ("area", "A(t): bubble area"),
        ("y_centroid", "yc(t): bubble centroid height"),
        ("D_eq", "D_eq(t): equivalent diameter"),
        ("wall_vapor_cells", "wall vapor cells"),
    ]

    for ax, (col, title) in zip(axes.ravel(), plots):
        ax.plot(df["step"], df[col], linewidth=2)
        ax.set_xlabel("step")
        ax.set_ylabel(col)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

        if departure_step is not None:
            ax.axvline(departure_step, linestyle="--", linewidth=1.2)

    fig.suptitle("Bubble metrics from phase field", fontsize=14)
    fig.tight_layout()

    path = os.path.join(OUTPUT_DIR, "bubble_metrics_4in1.png")
    plt.savefig(path, dpi=220)
    plt.close()

    print(f"保存：{path}")


# ============================================================
# 主程序
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = find_csv_files()

    rows = []
    for f in files:
        rows.append(compute_metrics_from_file(f))

    metrics = pd.DataFrame(rows)
    metrics = metrics.sort_values("step").reset_index(drop=True)

    out_csv = os.path.join(OUTPUT_DIR, "bubble_metrics.csv")
    metrics.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"保存：{out_csv}")

    departure_step = estimate_departure_step(metrics)

    if departure_step is not None:
        print(f"估计脱落时刻：step = {departure_step}")
    else:
        print("未能自动估计脱落时刻，可以根据动画人工判断。")

    save_single_plot(
        metrics,
        "step",
        "area",
        "bubble area A",
        "Bubble area A(t)",
        "bubble_area.png",
        departure_step
    )

    save_single_plot(
        metrics,
        "step",
        "y_centroid",
        "bubble centroid height yc",
        "Bubble centroid height yc(t)",
        "bubble_centroid_y.png",
        departure_step
    )

    save_single_plot(
        metrics,
        "step",
        "D_eq",
        "equivalent diameter D_eq",
        "Bubble equivalent diameter D_eq(t)",
        "bubble_equivalent_diameter.png",
        departure_step
    )

    save_single_plot(
        metrics,
        "step",
        "wall_vapor_cells",
        "wall vapor cells",
        "Wall vapor cells",
        "wall_vapor_cells.png",
        departure_step
    )

    save_4in1_plot(metrics, departure_step)

    print("全部分析完成。")


if __name__ == "__main__":
    main()