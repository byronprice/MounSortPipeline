#!/bin/bash

source activate mountainsort

for j in $( ls -d */); do
	cd $j
	cp ~/Documents/Current-Projects/MounSortPipeline/params.json .
	counter=0
	for i in $( ls raw.nt* ); do
		let counter=counter+1
		ml-run-process ephys.whiten --inputs timeseries:$i --outputs timeseries_out:pre.mda$counter.prv
		ml-run-process ms4alg.sort --inputs timeseries:pre.mda$counter.prv --outputs firings_out:firings$counter.mda --parameters detect_sign:-1 adjacency_radius:-1 detect_threshold:1 clip_size:40
        done
        cd ..
done

source deactivate mountainsort

# export PATH="/home/byron/MATLAB/R2017b/bin:$PATH"
# matlab -nodisplay -r "files = dir('*');for ii=1:length(files);if files(ii).isdir;cd(files(ii).name));firingsmda2mat(files(ii).name);cd ..;end;end;MounSortCompileData;exit"
