# Repeated measures ANOVA
I always forget how to do this and lose old scripts in the nether, so saving here.

In ```R```, read in your data, example below.
```
library(tidyverse)
library(emmeans)
library(nlme)

BEAKER_PC = read_csv("BEAKER_pc.csv") %>%
  mutate(across(c(group, Rep, Sample, Day), as.factor)) %>%
  mutate(PC_mgL = sqrt(PC_mgL)) %>%
  mutate(APC_mgL = sqrt(APC_mgL))

head(BEAKER_PC)

# A tibble: 6 × 7
  Sample Day   PC_mgL APC_mgL group    Rep   
  <fct>  <fct>  <dbl>   <dbl> <fct>    <fct>      
1 B4     0      0.515    9.05 BG11half BG11half_B4
2 B8     0      0.580    12.4 BG11o    BG11o_B8   
3 B7     0      0.598    2.60 BG11o    BG11o_B7   
4 B9     0      0.612    6.99 BG11o    BG11o_B9   
5 B5     0      0.442    5.69 BG11half BG11half_B5
6 B6     0      0.423    1.30 BG11half BG11half_B6
```
Confirm normality
```
ggpubr::ggqqplot(BEAKER_PC$PC_mgL)
ggpubr::ggqqplot(BEAKER_PC$APC_mgL)
```
Run the model 
```
model_BEAKER_PC <- lme(
  PC_mgL ~ Day*group,
  random = ~1 | Rep,
  data = BEAKER_PC
)
```
Do the pairwise comparisions, but select only for comparison of 'Day', not global comparision.
```
pairs(emmeans(model_BEAKER_PC, ~ group | Day))

Day = 0:
 contrast         estimate  SE df t.ratio p.value
 BG11 - BG11half     -3.40 189  6  -0.018  0.9998
 BG11 - BG11o        -5.38 189  6  -0.028  0.9996
 BG11half - BG11o    -1.99 169  6  -0.012  0.9999

Day = 7:
 contrast         estimate  SE df t.ratio p.value
 BG11 - BG11half     -6.11 169  6  -0.036  0.9993
 BG11 - BG11o        12.61 169  6   0.075  0.9969
 BG11half - BG11o    18.72 169  6   0.111  0.9933

Day = 14:
 contrast         estimate  SE df t.ratio p.value
 BG11 - BG11half    -96.21 169  6  -0.569  0.8408
 BG11 - BG11o       369.07 169  6   2.184  0.1527
 BG11half - BG11o   465.28 169  6   2.754  0.0740

Day = 28:
 contrast         estimate  SE df t.ratio p.value
 BG11 - BG11half    -90.92 169  6  -0.538  0.8561
 BG11 - BG11o       454.83 169  6   2.692  0.0800
 BG11half - BG11o   545.75 169  6   3.230  0.0410

Degrees-of-freedom method: containment 
P value adjustment: tukey method for comparing a family of 3 estimates 
```

