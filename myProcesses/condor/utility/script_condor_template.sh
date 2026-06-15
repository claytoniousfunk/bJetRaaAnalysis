#!/bin/bash

ENV_SETUP

PRERUN

EXECUTABLE

echo 'run finished'

OUTPUT

if [ $? -ne 0 ]
then
    echo 'ERROR: problem copying job directory back'
else
    echo 'job directory copy succeeded'
fi
