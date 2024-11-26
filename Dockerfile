# Tested on alpine 3.20
FROM alpine:latest AS builder
# install dependencies
RUN apk update \
    && apk upgrade \
    && apk add --no-cache \
    clang \
    clang-dev \
    make \
    cmake \
    boost-dev \
    hdf5-dev


RUN mkdir -p /app/src
# copy files
COPY CMakeLists.txt /app
COPY src /app/src

WORKDIR /app/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make -j 8

FROM alpine:latest
# Copy libraries from build stage
# COPY --from=builder /usr/lib/ /usr/lib/
# this copies much more than needed. It's easier to just reinstall only what is needed vs copying.
RUN apk update && \
    apk add --no-cache \ 
    libstdc++ \
    hdf5-dev
# Copy application from build stage
COPY --from=builder /app/build/nazcaBPM.x /app/nazcaBPM.x

ENTRYPOINT [ "/app/nazcaBPM.x" ]