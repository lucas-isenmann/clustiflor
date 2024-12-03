import numpy as np
import matplotlib.pyplot as plt

# Read data
data = []
for line in open('lol1.csv', 'r'):
    x, y = map(float, line.strip().split())
    data.append([x, y])

# Sort by x
data.sort(key=lambda x: x[0])
print(data)

# Create line graph
length = 0.5
h = 0.01

x_avg = []
y_avg = []

for a in range(int(length/h)):
    c = a*h
    s = 0
    n = 0
    for [x,y] in data:
        if x > c+h:
            break
        elif x > c-h:
            n += 1
            s += y
    if n > 0:
        mean = s/n
        print(c, mean)
        x_avg.append(c)
        y_avg.append(mean)

# x_avg = np.arange(min(data[:, 0]), max(data[:, 0]) + h, h)
# y_avg = np.array([np.mean([y for x, y in data if x >= xi - h and x <= xi + h])
#                   for xi in x_avg])

plt.figure(figsize=(10, 6))
plt.plot(x_avg, y_avg)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Averaged Line Graph')
plt.grid(True)
plt.show()
