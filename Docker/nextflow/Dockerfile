FROM openjdk:8-jre-alpine

RUN apk update && apk add bash coreutils curl graphviz python3 py3-pip

RUN pip3 install --upgrade pip
RUN pip3 install awscli

# see https://blogs.oracle.com/java-platform-group/java-se-support-for-docker-cpu-and-memory-limits
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow AWS_DEFAULT_REGION=us-east-1

# copy docker client
RUN mkdir /data/
COPY main.nf /data/
COPY ufl-apply_panels.nf /data/
COPY ufl-germline.nf /data/
COPY ufl-multiqc.nf /data/
COPY nextflow.config /data/
COPY bin/ /data/bin/
COPY modules/ /data/modules/
COPY dist/docker /usr/local/bin/docker
COPY entry.sh /usr/local/bin/entry.sh
COPY nextflow /usr/local/bin/nextflow

# download runtime
RUN mkdir /.nextflow \
&& touch /.nextflow/dockerized \
&& chmod 755 /usr/local/bin/entry.sh \
&& chmod 755 /usr/local/bin/nextflow \
&& nextflow info

# define the entry point
ENTRYPOINT ["/usr/local/bin/entry.sh"]
