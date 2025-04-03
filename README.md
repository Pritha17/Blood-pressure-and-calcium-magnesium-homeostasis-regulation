# Blood-pressure-and-calcium-magnesium-homeostasis-regulation
# About
This is a mathematical model of blood pressure regulation and calcium-magnesium homeostasis in a male rat. The model was used to simulate (i) hypertension by different hypertensive stimuli, (ii) dietary magnesium deficiency, dietary calcium deficiency, and vitamin D3 deficiency, and (iii) primary hyperparathyroidism.

# Instructions
1. The magnesium_mod.m file contains all the equations for calcium-magnesium homeostasis regulation and the params_male.txt file contains all the corresponding parameter values.
2. The Ca_Mg_bp_reg_mod.m file contains ll the equations for blood pressure regulation and the get_pars.m file contains all the corresponding parameter values.
3. To simuate different types of hypertension and plot the results, run sim_HTN.m file.
4. To simuate dietary magnesium deficiency, dietary calcium deficiency, and vitamin D3 deficiency and plot the results, run sim_Ca_Mg_D3_deficit.m file.
5. To simuate primary hyperparathyroidism and plot the results, run sim_HPTH.m file.
6. The Rat_Data folder contains the baseline and simulation datafiles.

# Related work
Please cite the foolowing paper when using this model.

* [2024 Dutta et al. "Modeling calcium and magnesium balance: Regulation by calciotropic hormones and adaptations under varying dietary intake"](https://www.cell.com/iscience/fulltext/S2589-0042(24)02302-2)
* [2020 Ahmed et al. "Sex-specific computational models for blood pressure regulation in the rat"](https://journals.physiology.org/doi/full/10.1152/ajprenal.00376.2019)
