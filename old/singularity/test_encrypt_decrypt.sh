#!/bin/bash

orig_key=$(mutein encrypt-key --debug --output-file tmp/test_key 2> tmp/orig_key)

decrypted_key=$(mutein decrypt-key --debug --input-file tmp/test_key --output-file tmp/test_out 2> tmp/dec_key)

if diff -q tmp/orig_key tmp/dec_key ; then
    echo "keys matched"
    exit 0
else
    echo "keys dont match!"
    exit 1
fi
