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

DATA_DIR = "AUTO"

NX = 160
NY = 100

DROP_TOP_ROW = True
MAX_FRAMES = None

FPS = 6
INTERVAL = 180

DRAW_INTERFACE = True
DRAW_QUIVER = True

QUIVER_STEP = 6
QUIVER_SCALE = 0.15

SAVE_SINGLE_GIFS = True
SAVE_COMBINED_GIF = True
SAVE_CROPPED_SPEED_GIF = True
SAVE_VORTICITY_GIF = True

# ------------------------------------------------------------
# 色标范围
# ------------------------------------------------------------

SPEED_VMIN = 0.0
SPEED_VMAX = 0.015

# 裁剪速度图，只看气泡活动区域，避免顶部速度带干扰
Y_MAX_SPEED_PLOT = 75

# 温度场色标
T_VMIN = 1.000
T_VMAX = 1.045

# rho 是由 phi 映射的可视化密度
RHO_VMIN = 0.12
RHO_VMAX = 1.00

# 涡量色标，None 表示自动用分位数
VORTICITY_ABS_MAX = None


# ============================================================
# 自动寻找 CSV 文件
# ============================================================

def find_data_files():
    if DATA_DIR == "AUTO":
        candidates = [
            os.path.join("csv_out", "fields_*.csv"),
            os.path.join(".", "fields_*.csv"),
        ]
    else:
        candidates = [
            os.path.join(DATA_DIR, "fields_*.csv")
        ]

    for pattern in candidates:
        found = glob.glob(pattern)
        if found:
            return pattern, found

    raise FileNotFoundError(
        "没有找到 fields_*.csv。\n"
        "请确认：\n"
        "1. CSV 文件名类似 fields_0000000.csv\n"
        "2. 如果脚本在工程根目录，CSV 应该在 csv_out 文件夹\n"
        "3. 如果脚本在 csv_out 里面，DATA_DIR 可以用 AUTO 或 '.'"
    )


def extract_step(filename):
    name = os.path.basename(filename)
    m = re.search(r"fields_(\d+)\.csv", name)
    if m:
        return int(m.group(1))
    return -1


pattern, files = find_data_files()
files = sorted(files, key=extract_step)

if MAX_FRAMES is not None:
    files = files[:MAX_FRAMES]

print(f"CSV 搜索路径：{pattern}")
print(f"找到 {len(files)} 个场文件")

if len(files) <= 1:
    print("警告：只找到 1 个 CSV 文件，生成的 GIF 只会有 1 帧。")
    print("请检查 csv_out 目录里是否有 fields_0000200.csv、fields_0000400.csv 等连续文件。")

print("前几个文件：")
for f in files[:5]:
    print(" ", os.path.basename(f))


# ============================================================
# 读取 CSV 为二维场
# ============================================================

def load_csv_as_arrays(csv_file):
    df = pd.read_csv(csv_file)

    required_cols = {"x", "y", "rho", "phi", "T", "dx", "dy", "speed", "solid"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{csv_file} 缺少列：{missing}")

    if DROP_TOP_ROW:
        df = df[df["y"] < NY - 1]

    ny_use = int(df["y"].max()) + 1

    arrays = {}

    for name in ["rho", "phi", "T", "dx", "dy", "speed", "solid"]:
        arr = np.full((ny_use, NX), np.nan, dtype=float)

        pivot = df.pivot(index="y", columns="x", values=name)

        y_index = pivot.index.to_numpy(dtype=int)
        x_index = pivot.columns.to_numpy(dtype=int)

        valid_y = (y_index >= 0) & (y_index < ny_use)
        valid_x = (x_index >= 0) & (x_index < NX)

        arr[np.ix_(y_index[valid_y], x_index[valid_x])] = pivot.to_numpy()[np.ix_(valid_y, valid_x)]

        arrays[name] = arr

    return arrays


# ============================================================
# 计算涡量
# omega = dv/dx - du/dy
# dx 是 u，dy 是 v
# ============================================================

def compute_vorticity(u, v):
    dvdx = np.zeros_like(v)
    dudy = np.zeros_like(u)

    dvdx[:, 1:-1] = 0.5 * (v[:, 2:] - v[:, :-2])
    dvdx[:, 0] = v[:, 1] - v[:, 0]
    dvdx[:, -1] = v[:, -1] - v[:, -2]

    dudy[1:-1, :] = 0.5 * (u[2:, :] - u[:-2, :])
    dudy[0, :] = u[1, :] - u[0, :]
    dudy[-1, :] = u[-1, :] - u[-2, :]

    return dvdx - dudy


# ============================================================
# 扫描色标范围
# ============================================================

def scan_range(field_name, percentile_clip=False):
    values = []

    for f in files:
        arr = load_csv_as_arrays(f)[field_name]
        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            values.append(arr)

    if not values:
        raise ValueError(f"无法扫描字段 {field_name} 的范围")

    values = np.concatenate(values)

    if percentile_clip:
        return np.percentile(values, 1), np.percentile(values, 99)

    return np.min(values), np.max(values)


def scan_vorticity_abs_max():
    vals = []

    for f in files:
        data = load_csv_as_arrays(f)
        omega = compute_vorticity(data["dx"], data["dy"])
        omega = omega[np.isfinite(omega)]
        if omega.size > 0:
            vals.append(omega)

    if not vals:
        return 0.01

    vals = np.concatenate(vals)
    return np.percentile(np.abs(vals), 99)


range_T_auto = scan_range("T")
range_rho_auto = scan_range("rho")
range_speed_auto = scan_range("speed", percentile_clip=True)

range_T = (
    range_T_auto[0] if T_VMIN is None else T_VMIN,
    range_T_auto[1] if T_VMAX is None else T_VMAX
)

range_rho = (
    RHO_VMIN if RHO_VMIN is not None else range_rho_auto[0],
    RHO_VMAX if RHO_VMAX is not None else range_rho_auto[1]
)

range_speed = (
    SPEED_VMIN,
    SPEED_VMAX
)

if VORTICITY_ABS_MAX is None:
    omega_abs_max = scan_vorticity_abs_max()
else:
    omega_abs_max = VORTICITY_ABS_MAX

if omega_abs_max <= 1.0e-12:
    omega_abs_max = 0.01

range_vorticity = (-omega_abs_max, omega_abs_max)

print("T auto range         =", range_T_auto)
print("rho auto range       =", range_rho_auto)
print("speed auto range     =", range_speed_auto)
print("T used range         =", range_T)
print("rho used range       =", range_rho)
print("speed used range     =", range_speed)
print("vorticity used range =", range_vorticity)


# ============================================================
# 删除旧 contour，兼容新版 Matplotlib
# ============================================================

def remove_contour(contour_obj):
    if contour_obj is None:
        return None

    try:
        contour_obj.remove()
    except Exception:
        try:
            for coll in contour_obj.collections:
                coll.remove()
        except Exception:
            pass

    return None


# ============================================================
# 单场动画：T / rho / speed
# ============================================================

def make_single_field_gif(field_name, gif_name, cmap, vmin=None, vmax=None, contour_color="black"):
    first = load_csv_as_arrays(files[0])
    field0 = first[field_name]
    phi0 = first["phi"]

    if vmin is None or vmax is None:
        vmin, vmax = scan_range(field_name)

    fig, ax = plt.subplots(figsize=(9, 5.8))

    img = ax.imshow(
        field0,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest"
    )

    cbar = plt.colorbar(img, ax=ax)
    cbar.set_label(field_name)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(0, NX - 1)
    ax.set_ylim(0, NY - 2 if DROP_TOP_ROW else NY - 1)

    title = ax.set_title("")
    contour_holder = {"obj": None}

    if DRAW_INTERFACE:
        try:
            contour_holder["obj"] = ax.contour(
                phi0,
                levels=[0.0],
                colors=contour_color,
                linewidths=1.0,
                origin="lower"
            )
        except Exception:
            contour_holder["obj"] = None

    def update(i):
        csv_file = files[i]
        step = extract_step(csv_file)
        data = load_csv_as_arrays(csv_file)

        img.set_data(data[field_name])

        if DRAW_INTERFACE:
            contour_holder["obj"] = remove_contour(contour_holder["obj"])
            try:
                contour_holder["obj"] = ax.contour(
                    data["phi"],
                    levels=[0.0],
                    colors=contour_color,
                    linewidths=1.0,
                    origin="lower"
                )
            except Exception:
                contour_holder["obj"] = None

        title.set_text(f"{field_name} field | step = {step}")
        return [img, title]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(files),
        interval=INTERVAL,
        blit=False,
        repeat=True
    )

    print(f"正在保存 {gif_name} ...")
    anim.save(gif_name, writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"{gif_name} 保存完成")


# ============================================================
# 裁剪版速度场动画
# 只看 y <= Y_MAX_SPEED_PLOT，避免顶部速度带影响
# ============================================================

def make_cropped_speed_gif(gif_name="velocity_speed_cropped.gif"):
    first = load_csv_as_arrays(files[0])

    y_max = min(Y_MAX_SPEED_PLOT, first["speed"].shape[0])

    speed0 = first["speed"][:y_max, :]
    phi0 = first["phi"][:y_max, :]

    fig, ax = plt.subplots(figsize=(9, 5.8))

    img = ax.imshow(
        speed0,
        origin="lower",
        cmap="plasma",
        vmin=range_speed[0],
        vmax=range_speed[1],
        aspect="auto",
        interpolation="nearest"
    )

    cbar = plt.colorbar(img, ax=ax)
    cbar.set_label("speed")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(0, NX - 1)
    ax.set_ylim(0, y_max - 1)

    title = ax.set_title("")
    contour_holder = {"obj": None}

    if DRAW_INTERFACE:
        try:
            contour_holder["obj"] = ax.contour(
                phi0,
                levels=[0.0],
                colors="white",
                linewidths=1.0,
                origin="lower"
            )
        except Exception:
            contour_holder["obj"] = None

    def update(i):
        csv_file = files[i]
        step = extract_step(csv_file)
        data = load_csv_as_arrays(csv_file)

        speed = data["speed"][:y_max, :]
        phi = data["phi"][:y_max, :]

        img.set_data(speed)

        if DRAW_INTERFACE:
            contour_holder["obj"] = remove_contour(contour_holder["obj"])
            try:
                contour_holder["obj"] = ax.contour(
                    phi,
                    levels=[0.0],
                    colors="white",
                    linewidths=1.0,
                    origin="lower"
                )
            except Exception:
                contour_holder["obj"] = None

        title.set_text(f"cropped speed field | step = {step}")
        return [img, title]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(files),
        interval=INTERVAL,
        blit=False,
        repeat=True
    )

    print(f"正在保存 {gif_name} ...")
    anim.save(gif_name, writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"{gif_name} 保存完成")


# ============================================================
# 涡量动画
# ============================================================

def make_vorticity_gif(gif_name="vorticity.gif"):
    first = load_csv_as_arrays(files[0])

    y_max = min(Y_MAX_SPEED_PLOT, first["speed"].shape[0])

    omega0 = compute_vorticity(first["dx"], first["dy"])[:y_max, :]
    phi0 = first["phi"][:y_max, :]

    fig, ax = plt.subplots(figsize=(9, 5.8))

    img = ax.imshow(
        omega0,
        origin="lower",
        cmap="seismic",
        vmin=range_vorticity[0],
        vmax=range_vorticity[1],
        aspect="auto",
        interpolation="nearest"
    )

    cbar = plt.colorbar(img, ax=ax)
    cbar.set_label("vorticity")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(0, NX - 1)
    ax.set_ylim(0, y_max - 1)

    title = ax.set_title("")
    contour_holder = {"obj": None}

    if DRAW_INTERFACE:
        try:
            contour_holder["obj"] = ax.contour(
                phi0,
                levels=[0.0],
                colors="black",
                linewidths=1.0,
                origin="lower"
            )
        except Exception:
            contour_holder["obj"] = None

    def update(i):
        csv_file = files[i]
        step = extract_step(csv_file)
        data = load_csv_as_arrays(csv_file)

        omega = compute_vorticity(data["dx"], data["dy"])[:y_max, :]
        phi = data["phi"][:y_max, :]

        img.set_data(omega)

        if DRAW_INTERFACE:
            contour_holder["obj"] = remove_contour(contour_holder["obj"])
            try:
                contour_holder["obj"] = ax.contour(
                    phi,
                    levels=[0.0],
                    colors="black",
                    linewidths=1.0,
                    origin="lower"
                )
            except Exception:
                contour_holder["obj"] = None

        title.set_text(f"vorticity field | step = {step}")
        return [img, title]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(files),
        interval=INTERVAL,
        blit=False,
        repeat=True
    )

    print(f"正在保存 {gif_name} ...")
    anim.save(gif_name, writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"{gif_name} 保存完成")


# ============================================================
# 综合动画：T + rho + speed/velocity
# ============================================================

def make_combined_gif(gif_name="combined_T_rho_flow.gif"):
    first = load_csv_as_arrays(files[0])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    axT, axRho, axU = axes

    imgT = axT.imshow(
        first["T"],
        origin="lower",
        cmap="inferno",
        vmin=range_T[0],
        vmax=range_T[1],
        aspect="auto",
        interpolation="nearest"
    )

    imgRho = axRho.imshow(
        first["rho"],
        origin="lower",
        cmap="viridis",
        vmin=range_rho[0],
        vmax=range_rho[1],
        aspect="auto",
        interpolation="nearest"
    )

    imgU = axU.imshow(
        first["speed"],
        origin="lower",
        cmap="plasma",
        vmin=range_speed[0],
        vmax=range_speed[1],
        aspect="auto",
        interpolation="nearest"
    )

    plt.colorbar(imgT, ax=axT, fraction=0.046, pad=0.04).set_label("T")
    plt.colorbar(imgRho, ax=axRho, fraction=0.046, pad=0.04).set_label("rho")
    plt.colorbar(imgU, ax=axU, fraction=0.046, pad=0.04).set_label("speed")

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(0, NX - 1)
        ax.set_ylim(0, NY - 2 if DROP_TOP_ROW else NY - 1)

    axT.set_title("Temperature T")
    axRho.set_title("Density rho")
    axU.set_title("Velocity magnitude")

    title = fig.suptitle("")

    contour_holders = [{"obj": None}, {"obj": None}, {"obj": None}]
    quiver_holder = {"obj": None}

    def draw_interfaces(data):
        if not DRAW_INTERFACE:
            return

        colors = ["cyan", "white", "white"]

        for ax, holder, color in zip(axes, contour_holders, colors):
            holder["obj"] = remove_contour(holder["obj"])
            try:
                holder["obj"] = ax.contour(
                    data["phi"],
                    levels=[0.0],
                    colors=color,
                    linewidths=1.0,
                    origin="lower"
                )
            except Exception:
                holder["obj"] = None

    def draw_quiver(data):
        if not DRAW_QUIVER:
            return

        if quiver_holder["obj"] is not None:
            try:
                quiver_holder["obj"].remove()
            except Exception:
                pass
            quiver_holder["obj"] = None

        yy, xx = np.mgrid[0:data["dx"].shape[0], 0:data["dx"].shape[1]]

        xs = xx[::QUIVER_STEP, ::QUIVER_STEP]
        ys = yy[::QUIVER_STEP, ::QUIVER_STEP]
        us = data["dx"][::QUIVER_STEP, ::QUIVER_STEP]
        vs = data["dy"][::QUIVER_STEP, ::QUIVER_STEP]

        quiver_holder["obj"] = axU.quiver(
            xs,
            ys,
            us,
            vs,
            color="white",
            scale=QUIVER_SCALE,
            scale_units="xy",
            angles="xy",
            width=0.0025
        )

    draw_interfaces(first)
    draw_quiver(first)

    def update(i):
        csv_file = files[i]
        step = extract_step(csv_file)
        data = load_csv_as_arrays(csv_file)

        imgT.set_data(data["T"])
        imgRho.set_data(data["rho"])
        imgU.set_data(data["speed"])

        draw_interfaces(data)
        draw_quiver(data)

        title.set_text(f"Pool boiling fields | step = {step}")
        return [imgT, imgRho, imgU, title]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(files),
        interval=INTERVAL,
        blit=False,
        repeat=True
    )

    print(f"正在保存 {gif_name} ...")
    anim.save(gif_name, writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"{gif_name} 保存完成")


# ============================================================
# 主程序
# ============================================================

if __name__ == "__main__":

    if SAVE_SINGLE_GIFS:
        make_single_field_gif(
            field_name="T",
            gif_name="temperature_T.gif",
            cmap="inferno",
            vmin=range_T[0],
            vmax=range_T[1],
            contour_color="black"
        )

        make_single_field_gif(
            field_name="rho",
            gif_name="density_rho.gif",
            cmap="viridis",
            vmin=range_rho[0],
            vmax=range_rho[1],
            contour_color="black"
        )

        make_single_field_gif(
            field_name="speed",
            gif_name="velocity_speed.gif",
            cmap="plasma",
            vmin=range_speed[0],
            vmax=range_speed[1],
            contour_color="white"
        )

    if SAVE_CROPPED_SPEED_GIF:
        make_cropped_speed_gif("velocity_speed_cropped.gif")

    if SAVE_VORTICITY_GIF:
        make_vorticity_gif("vorticity.gif")

    if SAVE_COMBINED_GIF:
        make_combined_gif("combined_T_rho_flow.gif")

    print("全部动画生成完成")