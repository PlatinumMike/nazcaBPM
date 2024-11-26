# TODO: switch to alpine build for smaller filesize. And/or a multi-stage build to save space.
FROM ubuntu:24.04
# install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends libhdf5-dev libboost-all-dev cmake make g++ && \
    apt-get clean

RUN mkdir -p /app/src
# copy files
COPY CMakeLists.txt /app
COPY src /app/src

WORKDIR /app/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make -j 8
ENTRYPOINT [ "./nazcaBPM.x" ]