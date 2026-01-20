# cherry-blossom-2026-internal-clock
Predicting 2026 peak bloom dates for five international locations using a transcriptomic proxy model (DAM4 gene activity) and Arrhenius thermal sum equations.

# 2026 International Cherry Blossom Prediction
**Name:** Bethany Gopinath  
**Methodology:** Internal Clock (Transcriptomic Proxy)

## Overview
This repository contains the 2026 peak bloom predictions for:
1. Washington D.C., USA
2. Kyoto, Japan
3. Liestal-Weideli, Switzerland
4. Vancouver, B.C., Canada
5. New York City, USA

## Methodology: The Internal Clock Model
Unlike linear GDD models, this project uses a two-phase physiological approach:
1. **Endodormancy (DAM4 Decay):** We simulate the activity of the $DAM4$ gene. We track hourly temperatures below 10.1°C; once 1,464 chilling hours (61 days) are reached, the tree is considered "awake."
2. **Ecodormancy (Arrhenius Forcing):** Post-wakeup, we calculate the Days of Thermal Sum (DTS) using the Arrhenius equation to predict the exact expansion of the buds.

## How to Reproduce
1. Install dependencies: `pip install -r requirements.txt`
2. Run the Quarto report: `quarto render reports/analysis.qmd`
3. View final results in `predictions/predictions.csv`.
