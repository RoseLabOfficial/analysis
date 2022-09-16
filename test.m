x = [5, 10, 20, 40, 60, 80];
y = [0, 1, 2, 4, 6, 4];
ymax = max(y, [], 2);
ydrop = ymax*0.707;
x1 = interp1(y, x, ydrop, "makima");