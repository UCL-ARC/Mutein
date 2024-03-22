#using the current Ubuntu LTS
FROM ubuntu:22.04

#use bash not sh
SHELL ["/bin/bash", "-c"]

#name of the limited user within the container
ARG CONTAINER_USER

#where under the users home folder to clone the pipeline repository to
ARG REPO_FOLDER

#where within the container will the external data folder of the TRE be mounted
#ARG GUEST_DATA_FOLDER

#url to download mamba install script from and name of the install script itself
ARG MAMBA_URL
ARG MAMBA_SH

#github host keys
# ARG github_ed25519="github.com ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIOMqqnkVzrm0SdG6UOoqKLsabgH5C9okWi0dh2l9GKJl"
# ARG github_ecdsa="github.com ecdsa-sha2-nistp256 AAAAE2VjZHNhLXNoYTItbmlzdHAyNTYAAAAIbmlzdHAyNTYAAABBBEmKSENjQEezOmxkZMy7opKgwFB9nkt5YRrYMjNuG5N87uRgg6CLrbo5wAdT/y6v0mKV0U2w0WZ2YB/++Tpockg="
# ARG github_rsa="github.com ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQCj7ndNxQowgcQnjshcLrqPEiiphnt+VTTvDP6mHBL9j1aNUkY4Ue1gvwnGLVlOhGeYrnZaMgRK6+PKCUXaDbC7qtbW8gIkhL7aGCsOr/C56SJMy/BCZfxd1nWzAOxSDPgVsmerOBYfNqltV9/hWCqBywINIR+5dIg6JTJ72pcEpEjcYgXkE2YEFXV1JHnsKgbLWNlhScqb2UmyRkQyytRLtL+38TGxkxCflmO+5Z8CSSNY7GidjMIZ7Q4zMjA2n1nGrlTDkzwDCsw+wqFPGQA179cnfGWOWRVruj16z6XyvxvjJwbz0wQZ75XK5tKSb7FNyeIEs4TT4jk+S4dhPeAUC5y+bDYirYgM4GC7uEnztnZyaVWQ7B381AK4Qdrwt51ZqExKbQpTUNn+EjqoTwvqNj4kqx5QUCI0ThS/YkOxJCXmPUWZbhjpCg56i+2aB6CmK2JGhn57K5mj0MNdBXA4/WnwH6XoPWJzK5Nyu2zB3nAZp+S5hpQs+p1vN1/wsjk="

#add limited user account with an .ssh folder
RUN useradd --create-home --shell /bin/bash $CONTAINER_USER
RUN mkdir /home/$CONTAINER_USER/.ssh
RUN chown $CONTAINER_USER:$CONTAINER_USER /home/$CONTAINER_USER/.ssh
RUN chmod go-rwx /home/$CONTAINER_USER/.ssh

#install git
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y git wget
RUN rm -rf /var/lib/apt/lists/*

#switch to operating as the limited user
USER $CONTAINER_USER
WORKDIR /home/$CONTAINER_USER

#add github's host keys (is this needed for an https clone?)
# RUN echo "$github_ed25519" >> ~/.ssh/known_hosts
# RUN echo "$github_ecdsa" >> ~/.ssh/known_hosts
# RUN echo "$github_rsa" >> ~/.ssh/known_hosts

#clone the public pipeline repo
RUN mkdir -p $REPO_FOLDER
WORKDIR /home/$CONTAINER_USER/$REPO_FOLDER
RUN git clone https://github.com/UCL-ARC/Mutein.git
WORKDIR /home/$CONTAINER_USER/$REPO_FOLDER/Mutein
RUN git checkout docker4tre

#install mamba
WORKDIR /home/$CONTAINER_USER
RUN mkdir -p mamba
WORKDIR /home/$CONTAINER_USER/mamba
RUN wget $MAMBA_URL/$MAMBA_SH
RUN chmod u+x $MAMBA_SH
RUN ./$MAMBA_SH -b -p /home/$CONTAINER_USER/mambaforge
ENV PATH=/home/treuser/mambaforge/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN conda init
RUN conda config --set auto_activate_base false

#install packages
WORKDIR /home/$CONTAINER_USER/$REPO_FOLDER/Mutein
RUN ./mutein/scripts/setup_conda_envs.sh mutein/conda_envs/*.yml

#change into same folder that run_container.sh launches from
#for convenience
WORKDIR /home/$CONTAINER_USER/$REPO_FOLDER/Mutein/docker

#set default command to bash
CMD ["bash"]
