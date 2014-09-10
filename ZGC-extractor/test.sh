#!/bin/bash

./rc < test.in > test.debug 2> test.out
if diff -bu test.out test.out.gold ; then
    echo PASSED
else
    echo FAILED
fi
