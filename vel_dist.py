import matplotlib.pyplot as plt

with open("al_after.cfg", "r") as f:
    line = [s.split() for s in f.readlines()]
    vel_list = [[float(a) for a in v[3:6]] for v in line if len(v) == 6]

print(len(vel_list))
f, ax = plt.subplots(4, 1)
for i, v_comp in enumerate(zip(*vel_list)):
    ax[i].hist(v_comp, bins=15)

ax[3].hist([sum([x**2 for x in v]) ** 0.5 for v in vel_list], bins=15)

f.show()
# f.savefig("vel_dist.png")

