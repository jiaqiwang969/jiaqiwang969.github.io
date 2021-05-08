#!/bin/bash
for file in `ls *.h`;do mv $file `echo $file|sed 's/_0_T//g'`;done;
