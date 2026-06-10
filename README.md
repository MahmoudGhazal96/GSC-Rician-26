# GSC-Rician

Code accompanying the paper "Rare-Event Simulation of Outage Probability in GSC/MRC Systems under Rician Fading" (Ghazal, Ben Rached, Al-Naffouri). It reproduces the outage probability results in Section VI.



## Files

1. `PIS_ET_CE.m` can be used to generate Table III. It contains the competitive MC methods for the case of independent fading and equal variance: partition importance sampling (PIS), exponential twisting (ET), and cross-entropy (CE).
2. `ConstantCorrelationMethods.m` can be used to generate Table VII. It contains the enhanced MC methods that work for the constant correlation (and equal means) case: cross-entropy correlated (CEC), and exponential twisting correlated (ETC).
3. `ETCvsSong.m` can be used to generate Table VIII. It contains the only method in the paper that works for arbitrary correlation structures: ETC. It is compared with an asymptotic expression derived by Song, and with naive Monte-Carlo (NMC) for viable probabilities.

## Note
MATLAB code in this repository was edited with assistance from Claude and ChatGPT for clarity and structure.
