#!/bin/bash

source activate mountainsort

cp ~/Documents/Current-Projects/MounSortPipeline/params.json .

ml-run-process ephys.bandpass_filter --inputs timeseries:raw.full.mda --outputs timeseries_out:filt.mda.prv --parameters samplerate:30000 freq_min:300 freq_max:6000
ml-run-process ephys.whiten --inputs timeseries:filt.mda.prv --outputs timeseries_out:pre.mda.prv
ml-run-process ms4alg.sort --inputs timeseries:pre.mda.prv --outputs firings_out:rawfire.mda --parameters detect_sign:-1 adjacency_radius:-1 detect_threshold:1

source deactivate mountainsort
