Structure selection in time-varying coefficient Cox model
=========================================================

In time-varying coefficient Cox model, it is of great practical interest to accurately identify the structure of covariate effects, covariates with null effect, constant effect and truly time-varying effect, and estimate the corresponding
regression coefficients. Combining the ideas of local polynomial smoothing and group nonnegative garrote, we develop a new penalization approach to achieve such goals. See the [paper for details](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987133/)

### What is inside
File **KGNG.R** contains the main function *KGNG* to calculate KGNG (kernel group nonnegative garrote estimator) or KGNG2 (KGNG with preliminary step) in time-varying coefficient Cox model. This method can automatically identify the structure of covariates, i.e., covariates with null effect(O), covariates with constant effect (C) and covariates with truly time-varying effect (NC), and estimate the corresponding regression coefficients. The structure of covariate effects has been automatically identified. 

File **simPBC.R** applies the *KGNG* method to analyze the PBC data (Fleming and Harrington (1991)). The results are shown in the following Figures and Table.

|![Image of convexlar ls](/figures/Figure3.png)|![Image of convexlasso ls](/figures/Figure4.png)|
|:---:|:---:|

![Image of convexlar ls](/figures/Table4.png)<!-- .element height="20%" width="20%" -->

About the PBC data: the data is from the Mayo Clinic trial in primary biliary cirrhosis (PBC) of the liver conducted between 1974 and 1984. The primary biliary cirrhosis is a chronic disease in which the bile ducts in one’s liver are slowly destroyed. In the study, 312 out of 424 patients who participated in the randomized trial were eligible for the analysis. There are 17 covariates: trtmt=treatment (Yes/No), age (in 10 years), gender=female/male, ascites=presence of ascites (Yes/No), hypato=presence of hepatomegaly (Yes/No), spiders=presence of spiders, edema=severity of oedema (0 denotes no oedema,
0.5 denotes untreated or successfully treated oedema and 1 denotes unsuccessfully treated oedema), logbili=logarithm of serum bilirubin (mg/dl), chol=serumcholesterol (mg/dl), logalb=logarithm of albumin (gm/dl), copper=urine copper (mg/day), alk=alkaline phosphatase (U/liter), sgot=liver enzyme (U/ml), trig=triglicerides (mg/dl), platelet=platelets per 10−3 ml3, logprotime=logarithm of prothrombin time (seconds), stage=histologic stage of disease (category: 1, 2, 3 or 4). 

### Citation
```
@article{xiao2016joint,
  title={Joint structure selection and estimation in the time-varying coefficient Cox model},
  author={Xiao, Wei and Lu, Wenbin and Zhang, Hao Helen},
  journal={Statistica Sinica},
  volume={26},
  number={2},
  pages={547},
  year={2016},
  publisher={NIH Public Access}
}
```
### Contacts
Wei Xiao <wxiao0421@gmail.com>         


