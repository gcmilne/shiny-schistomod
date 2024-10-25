# Shiny app for simulating a schistosomiasis transmission model

## Published scientific paper:
Gregory C. Milne, Rebecca C. Oettle, Charles Whittaker, Narcis B. Kabaterine, Maria-Gloria Basáñez, Joanne P. Webster, Martin Walker, Shona Wilson, 2024. Revisiting immunity versus exposure in schistosomiasis: a mathematical modelling study of delayed concomitant immunity. *PNAS Nexus*, https://doi.org/10.1093/pnasnexus/pgae471.

## User guide:
The app is run using a publicly available model, which can be installed as an R package found [here](https://github.com/gcmilne/SchistoTransmissionModel).

The user can simulate different mass drug adminstration (MDA) regimens by varying:
- immuno-epidemiological parameters (e.g. R0 & the strength of acquired immunity)
- the total number of rounds of MDA
- the frequency of MDA
- the MDA coverage in school-aged children (5-15 year olds)
- the MDA coverage in adults (16+ year olds)

The app outputs time-specific means of the following model outputs (which can be toggled on or off via the user interface):
- epg: the eggs per gram of faeces, a measure of the intensity of infection
- worm burden: the worm burden that underlies this infection intensity (i.e., the mean number of adult worms harboured by the population)
- prevalence: the proportion of the population with eggs in faeces
- smTAL1-IgE: the optical density of IgE antibodies raised against a Schistosoma mansoni tegument protein which is released upon worm death (both naturally and through MDA)

The user can also specify for how many years after the final MDA event the model output(s) should be shown, with options ranging from 5-20 years.
