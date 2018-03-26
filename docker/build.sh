#!/bin/bash

VERSION=0.0.2

sed -i '' "s/version \".*\"/version \"$VERSION\"/" Dockerfile

docker build --rm -t dongli/iap-cgfd-adv-cases:$VERSION . && \
docker push dongli/iap-cgfd-adv-cases:$VERSION
