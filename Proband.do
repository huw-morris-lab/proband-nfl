replace Nfl_meanconc=500 if Nfl_meanconc>500
drop if change_diag_ind
generate Nfl_meanconc=log( Nfl_meanconc)
logit change age_V1 gender Nfl_meanconc
lroc, nograph
recode hoehnyahrlast 0/3=0 3/max=1, generate(HighHY)
recode MOCAlast 0/17=1 18/max=0, generate(MOCAlastlow)
egen oddgroup = group( HighHY MOCAlastlow Death )
recode oddgroup 1=0 2=1 3=1 4=1 5=1 6=1 7=1
logit oddgroup age_V1 gender UPDRS80Axial, nolog
lroc, nograph
predict xb1, xb
logit oddgroup age_V1 gender logNfl , nolog
lroc, nograph
predict xb2, xb
roccomp oddgroup xb1 xb3
stset Survival_time, failure(oddgroup==1)
stcox gender age_V1 UPDRS80Axial
estat concordance
generate ranord=runiform()
generate long oldord=_n
sort ranord, stable
generate testset=mod(_n,2)
sort oldord
stcox gender age_V1 UPDRS80Axial
predict hr1
generate invhr1=1/hr1
stcox gender age_V1 logNfl
predict hr2
generate invhr2=1/hr2
generate censind=1-_d if _st==1
somersd _t invhr1 invhr2 invhr3 if _st==1 & testset==1, cenind(censind) tdist transf(c)
lincom invhr1-invhr2
reshape long MOCA hoehnyahr UPDRStotal UPDRSIII Schwab PDQ dyskinesias semanflu diseaseduration , i( Name ) j (year)
gen lnNfl_meanconcXdd = lnNfl_meanconc* diseaseduration
mixed MOCA diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
mixed hoehnyahr diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
mixed UPDRSIII diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
mixed Schwab diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
mixed PDQ diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
mixed semanflu diseaseduration lnNfl_meanconc lnNfl_meanconcXdd || Name: diseaseduration
logit oddgroup age_V1 gender, nolog
lroc, nograph
logit oddgroup age_V1 gender APOE_3 , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 , nolog
lroc, nograph
logit oddgroup age_V1 gender GBA  , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 UPDRS80Axial V1_seman_flu_score , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score Nfl_meanconc , nolog
logit oddgroup age_V1 gender rs429358 UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
drop if V1_hoehn_and_yahr >3
drop if V1_MOCA_total_adj <17
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score Nfl_meanconc , nolog
logit oddgroup age_V1 gender rs429358 UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial V1_seman_flu_score, nolog
lroc, nograph
drop if prob_of_pd_last <90
drop if change_diag_ind
drop if V1_hoehn_and_yahr >=3
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score Nfl_meanconc , nolog
logit oddgroup age_V1 gender rs429358 UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial, nolog
lroc, nograph
logit oddgroup age_V1 gender V1_seman_flu_score , nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial V1_seman_flu_score, nolog
lroc, nograph
logit oddgroup age_V1 gender UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender rs7412 UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender GBA UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
logit oddgroup age_V1 gender rs429358 GBA UPDRS80Axial V1_seman_flu_score std_lnNfl_meanconc , nolog
lroc, nograph
senspec Group Nfl_meanconc , sensitivity(sens) specificity(spec)
gen youden= sens-(1-spec)
egen youdenmax= max(youden), by(gender)
gen dist = sqrt((1-sens)^2 + (1-spec)^2)
egen distmin = min(dist), by(gender)
list sens spec  youdenmax   dist Nfl_meanconc if abs(youden -youdenmax)<0.0001
list sens spec  youden     distmin  Nfl_meanconc if abs(dist - distmin)<0.0001
xtile QNFL= Nfl_meanconc, nq(4)
stcox Nfl_meanconc, strata(QNFL)
stcox QNFL


