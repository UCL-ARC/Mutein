#using the current Ubuntu LTS
FROM ubuntu:22.04

#name of the limited user within the container
ARG username=treuser

#where under the users home folder to clone the pipeline repository to
ARG repo_path=repos

#where within the container will the external data folder of the TRE be mounted
ARG data_path=/shared/mutein_data

#add limited user account with .ssh folder
RUN useradd --create-home --shell /bin/bash $username
RUN mkdir /home/$username/.ssh
RUN chown $username:$username /home/$username/.ssh
RUN chmod go-rwx /home/$username/.ssh

#add github access ssh key
ADD secrets/id_ed25519 /home/$username/.ssh/id_ed25519
ADD secrets/id_ed25519.pub /home/$username/.ssh/id_ed25519.pub

RUN chmod go-wrx /home/$username/.ssh/id_ed25519
RUN chmod go-wx /home/$username/.ssh/id_ed25519.pub

USER $username
WORKDIR /home/$username

CMD ["bash"]
