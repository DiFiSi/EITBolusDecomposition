# EITBolusDecomposition
Demo code for gamma-variate-based non-negative matrix factorization (NMF) of electrical impedance tomography (EIT) bolus recordings used in our paper [Contrast-enhanced EIT robustly tracks regional lung perfusion compared to non-enhanced EIT and pulmonary CT](https://ieeexplore.ieee.org/abstract/document/10933520).

## Algorithm
Example (real) data is provided in ./data.mat, and the main script to showcase the algorithm is ./go.m. The script showcases:
1. Pre-processing algorithm - drift removal using the Heron onset detection algorithm and pre-bolus detrending
2. NMF with N = 3 (three compartments) - decomposes the entire signal into pre-lung, lung, and post-lung compartments on the basis of their gamma-variate-parametrized time dynamics
3. Lung splitting - split spatial map of lung compartment into right and left lung compartment spatial maps
4. NMF with N = 4 (four compartments) - decomposes the entire signal into right/left heart/lung compartments using spatial priors obtained in 3.

Overview of the algorithm (taken from paper)
![image](https://github.com/user-attachments/assets/814f4bbd-45a6-47f8-9494-897d38d9c6ee)

## Files
* NMF and gamma-variate model functions are found in ./model/
* Pre-processing functions are found in ./preprocess/
* Useful and third party functions are found in ./utils/ and ./others/ respectively 
