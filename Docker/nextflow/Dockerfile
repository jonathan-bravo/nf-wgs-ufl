FROM openjdk:8-jre-alpine

RUN apk update && apk add bash coreutils curl graphviz python3 py3-pip

RUN pip3 install --upgrade pip
RUN pip3 install awscli

# see https://blogs.oracle.com/java-platform-group/java-se-support-for-docker-cpu-and-memory-limits
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow AWS_DEFAULT_REGION=us-east-1

# copy docker client
COPY main.nf ~/
COPY ufl-apply_panels.nf ~/
COPY ufl-germline.nf ~/
COPY ufl-multiqc.nf ~/
COPY nextflow.config ~/
COPY bin/ ~/bin/
COPY modules/ ~/modules/
COPY dist/docker /usr/local/bin/docker
COPY entry.sh /usr/local/bin/entry.sh
COPY nextflow /usr/local/bin/nextflow

# download runtime
RUN mkdir /.nextflow \
 && touch /.nextflow/dockerized \
 && chmod 755 /usr/local/bin/nextflow \
 && chmod 755 /usr/local/bin/entry.sh \
 && nextflow info

# define the entry point
ENTRYPOINT ["/usr/local/bin/entry.sh"]
