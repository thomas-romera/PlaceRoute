import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("coloquinte_timing.csv")
detail = df[df.sub_step != "total"]

fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=False)

# --- UB sub-steps ---
df_ub = detail[detail.phase == "UB"].pivot_table(
    index="step", columns="sub_step", values="time_ms", aggfunc="sum"
)
df_ub.plot.bar(stacked=True, ax=axes[0], title="runUB sub-steps per iteration")
axes[0].set_ylabel("time (ms)")
axes[0].set_xlabel("")

# --- LB sub-steps ---
df_lb = detail[detail.phase == "LB"].pivot_table(
    index="step", columns="sub_step", values="time_ms", aggfunc="sum"
)
df_lb.plot.bar(stacked=True, ax=axes[1], title="runLB sub-steps per iteration")
axes[1].set_ylabel("time (ms)")
axes[1].set_xlabel("")

# --- Total per iteration (UB + LB totals stacked) ---
df_tot = df[df.sub_step == "total"].pivot_table(
    index="step", columns="phase", values="time_ms", aggfunc="sum"
)
df_tot.plot.bar(stacked=True, ax=axes[2], title="Total time per iteration (UB + LB)")
axes[2].set_ylabel("time (ms)")
axes[2].set_xlabel("step")

plt.tight_layout()
plt.savefig("timing.png", dpi=150)
plt.show()
