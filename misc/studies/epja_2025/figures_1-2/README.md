# 1-body local density plots for 16O and 240Pu deformed solutions

1. generate the HFB solutions with the commands

```bash
./hfb3 16O_deformed.hfb3
./hfb3 240Pu_deformed.hfb3
```

2. source a Python venv if needed
3. generate the .eps figures with the command

```bash
./plotLocalDensity.py
```
