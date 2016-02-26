#!/bin/bash

for file in ./out/*.bmp
do
	echo $file
	convert $file ./out/`basename $file .bmp`.png
	rm -f $file
done

mencoder mf://./out/* -mf fps=15 -ovc lavc -lavcopts vcodec=mpeg4 -o out.mpeg;

rm -f ./out/*
