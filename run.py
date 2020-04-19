import os
from pathlib import Path

assignment = 2
program = "rt168"

os.system("cmake .. && make -j4")
os.system("rm -f *.png *.zip")

pathlist = Path(f"../scenes/hw{assignment}").glob("*.test")

for path in pathlist:
    os.system(f"./{program} {path}")

os.system(f"zip homework{assignment}.zip *.png")