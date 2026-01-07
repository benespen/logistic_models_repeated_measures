# logistic_models_repeated_measures

Combining Bayesian logistic regression with mixed models for repeated measures

The mmrm package points to Mallinckrodt, Lane, Schnell, Peng and Mancuso (2008) 
<doi:10.1177/009286150804200402> for a review of multi-level or mixed models
in the context of clinical trials.

Mallinckrodt is using continuous endpoints, but it is also common to have binary
yes/no endpoints for clinical trials, which necessitates logistic regression
or another similar model.

Statistical Rethinking 2nd ed. by Richard McElreath offers examples of 
Bayesian logistic regression in chapters 13 and 14. This repository is an 
attempt to think through using logistic regression on repeated measures.
