#using the current Ubuntu LTS
FROM ubuntu:22.04

#name of the limited user within the container
ARG username=treuser

#where under the users home folder to clone the pipeline repository to
ARG repo_path=repos

#where within the container will the external data folder of the TRE be mounted
ARG data_path=/shared/mutein_data

#github host keys
ARG github_ed25519="github.com ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIOMqqnkVzrm0SdG6UOoqKLsabgH5C9okWi0dh2l9GKJl"
ARG github_ecdsa="github.com ecdsa-sha2-nistp256 AAAAE2VjZHNhLXNoYTItbmlzdHAyNTYAAAAIbmlzdHAyNTYAAABBBEmKSENjQEezOmxkZMy7opKgwFB9nkt5YRrYMjNuG5N87uRgg6CLrbo5wAdT/y6v0mKV0U2w0WZ2YB/++Tpockg="
ARG github_rsa="github.com ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQCj7ndNxQowgcQnjshcLrqPEiiphnt+VTTvDP6mHBL9j1aNUkY4Ue1gvwnGLVlOhGeYrnZaMgRK6+PKCUXaDbC7qtbW8gIkhL7aGCsOr/C56SJMy/BCZfxd1nWzAOxSDPgVsmerOBYfNqltV9/hWCqBywINIR+5dIg6JTJ72pcEpEjcYgXkE2YEFXV1JHnsKgbLWNlhScqb2UmyRkQyytRLtL+38TGxkxCflmO+5Z8CSSNY7GidjMIZ7Q4zMjA2n1nGrlTDkzwDCsw+wqFPGQA179cnfGWOWRVruj16z6XyvxvjJwbz0wQZ75XK5tKSb7FNyeIEs4TT4jk+S4dhPeAUC5y+bDYirYgM4GC7uEnztnZyaVWQ7B381AK4Qdrwt51ZqExKbQpTUNn+EjqoTwvqNj4kqx5QUCI0ThS/YkOxJCXmPUWZbhjpCg56i+2aB6CmK2JGhn57K5mj0MNdBXA4/WnwH6XoPWJzK5Nyu2zB3nAZp+S5hpQs+p1vN1/wsjk="

#add limited user account with an .ssh folder
RUN useradd --create-home --shell /bin/bash $username
RUN mkdir /home/$username/.ssh
RUN chown $username:$username /home/$username/.ssh
RUN chmod go-rwx /home/$username/.ssh

#install git
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y git
RUN rm -rf /var/lib/apt/lists/*

#switch to operating as the limited user
USER $username
WORKDIR /home/$username

#add github's host keys
RUN echo "$github_ed25519" >> ~/.ssh/known_hosts
RUN echo "$github_ecdsa" >> ~/.ssh/known_hosts
RUN echo "$github_rsa" >> ~/.ssh/known_hosts

#clone the public pipeline repo into the container
RUN mkdir -p $repo_path
WORKDIR /home/$username/$repo_path
RUN git clone https://github.com/UCL-ARC/Mutein.git
WORKDIR /home/$username/$repo_path/Mutein
RUN git checkout robertvi_tre_docker

CMD ["bash"]
