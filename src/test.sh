#!/bin/sh

mkdir test
cd test
file_name="output.txt"
julia ../hello.jl > $file_name 
echo "I wrote a file called $file_name"