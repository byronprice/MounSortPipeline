#!/bin/bash

for j in $( ls -d */); do
	cd $j
	cp ~/Documents/Current-Projects/MounSortPipeline/mountainsort3.mlp .
	cp ~/Documents/Current-Projects/MounSortPipeline/params.json .
	counter=0
	for i in $( ls raw.nt* ); do
		let counter=counter+1
		mlp-run mountainsort3.mlp sort --raw=$i --firings_out=firings$counter.mda --_params=params.json --curate=true
        done
        cd ..
done

export PATH="/home/byron/MATLAB/R2017b/bin:$PATH"
matlab -nodisplay -r "files = dir('*');for ii=1:length(files);if files(ii).isdir;cd(files(ii).name));firingsmda2mat(files(ii).name);cd ..;end;end;MounSortCompileData;exit"
