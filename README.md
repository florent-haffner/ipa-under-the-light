# IPA Under the light - XAI

## Getting started
Information : y is already normalized (mean 0 & std 1), the scaler will not be given.

### Python
Follow the next commands:
```
conda create -n ipa-under-the-light python=3.10

conda activate ipa-under-the-light

conda install pytorch=2.1 pytorch-cuda=11.8 -c pytorch -c nvidia

pip install -r requirements.txt
```

### Matlab
Install [MATLAB R2022b](https://www.mathworks.com/products/new_products/release2022b.html) and [Eigenvector PLS_Toolbox version 9.2](https://wiki.eigenvector.com/index.php?title=Release_Notes_Version_9_2).

## Using the software

### Python

```bash
conda activate ipa-under-the-light

cd notebooks

jupyter lab .
```

### Matlab
Open the script within the eponymous software.
