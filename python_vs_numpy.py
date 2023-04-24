import numpy as np
from timeit import timeit

size = 100000000

def fill_python_array(arr):
    for i in range(size):
        arr[i] = "C"

def fill_numpy_array(arr):
    for i in range(size):
        arr[i] = 123

python_durations = []
numpy_durations = []

for i in range(0, 7):
    print("Runs", i+1, "of", 7)
    python_bytes = ["A"] * size
    numpy_bytes = np.zeros(size, dtype=np.uint8)
    duration = timeit(lambda: fill_python_array(python_bytes), number=1)
    python_durations.append(round(duration * 1000))
    duration = timeit(lambda: fill_numpy_array(numpy_bytes), number=1)
    numpy_durations.append(round(duration * 1000))

print("Python:", sorted(python_durations))
print("Numpy: ", sorted(numpy_durations))
