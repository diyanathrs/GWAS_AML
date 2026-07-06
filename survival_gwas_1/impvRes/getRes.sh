#!/bin/bash

while read link;
do
    echo -e "${link}"
    wget ${link}

done < reslink
