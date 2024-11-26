# TODO: Use a multi-stage build to save space.
# Tested on alpine 3.20
FROM alpine:latest
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
ENTRYPOINT [ "./nazcaBPM.x" ]