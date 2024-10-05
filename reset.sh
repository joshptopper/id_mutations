#!/bin/bash

rm *amb
rm *ann
rm *bwt
rm *pac
rm *sa
rm *fai
cd fastqs
rm *trim*
cd ..
rmdir fastqs
cd bams
rm *sort*
cd ..
rmdir bams
rm report*