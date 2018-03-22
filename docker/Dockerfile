FROM centos:7

LABEL description="This is an environment for running advection cases."
LABEL maintainer="Li Dong <dongli@lasg.iap.ac.cn>"
LABEL version="0.0.1"

RUN yum -y update && \
    yum -y install epel-release https://centos7.iuscommunity.org/ius-release.rpm && \
    yum -y install gcc-gfortran gsl-dev cmake3 make git2u wget && \
    yum -y install netcdf-fortran-devel && \
    yum clean all && rm -rf /var/cache/yum

RUN yum -y install libX11-devel cairo-devel && \
    cd /tmp && wget http://static.longrunweather.com/packages/ncl_ncarg-6.4.0-CentOS7.3_64bit_nodap_gnu485.tar.gz && \
    cd /usr/local && tar xf /tmp/ncl_ncarg-6.4.0-CentOS7.3_64bit_nodap_gnu485.tar.gz && \
    rm /tmp/ncl_ncarg-6.4.0-CentOS7.3_64bit_nodap_gnu485.tar.gz && \
    yum -y erase libX11-devel cairo-devel && \
    yum clean all && rm -rf /var/cache/yum

RUN cd /usr/include && \
    ln -s /usr/lib64/gfortran/modules/netcdf.mod . && \
    ln -s /usr/lib64/gfortran/modules/typesizes.mod . && \
    cd /home && \
    git clone https://github.com/dongli/iap-cgfd-adv-cases && \
    cd iap-cgfd-adv-cases && \
    mkdir build && \
    cd build && \
    cmake3 .. && \
    make

WORKDIR /home/iap-cgfd-adv-cases/build
ENV NCARG_ROOT=/usr/local
