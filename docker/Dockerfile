FROM ubuntu:22.04

ARG username=treuser

RUN \
    useradd --create-home --shell /bin/bash $username

USER $username
WORKDIR /home/$username

CMD ["bash"]
