import matplotlib.pyplot as plt

with open("strain.txt", "r") as f:

    str_eng = [[float(v) for v in line.split(",")[1:3]] for line in f.readlines()]

slope = []
for i in range(1, len(str_eng) - 1):
    strl, engl = str_eng[i - 1]
    strr, engr = str_eng[i + 1]
    slp = (engr - engl)/(strr - strl)
    if slp < 0: 
        continue
    slope.append((str_eng[i][0], slp))


plt.plot(*zip(*slope), marker="+", markersize=6, linestyle="")
plt.savefig("s-s.png")



