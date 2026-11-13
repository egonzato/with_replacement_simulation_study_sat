# with_replacement_simulation_study_sat

when few patients are available in the historical control pool some treated participants may not find a match. throwing away patients translates to throwing away important information and possibly making the wrong decision in terms of testing. can we fix this problem if performing matching with replacement?

# done

- dgm is based on previous simulation study, see `simulation_study_sat`
- adapted matching functions to allow for `replacement=TRUE`
- add weights
    - number of controls in cluster should sum to one (`cweight`)
    - number of times a unique control is used should sum to one (`iweight`)

# to do

- write shuffling function (both for RR and HR)
- write variance estimator function (Austin et al. adapted from Miglioretti et al.)
- think about possible solutions
    - adapt variance estimator to use weights
    - jacknife weights
    - make functions use single controls not more than in z% of clusters
- decide which estimand to go for (main issue arised is variance estimation)
    - hazard ratio not really meaningful (used in Austin et al.)
    - risk ratio needs bootstrap (but bootstrap is not valid for this case)
