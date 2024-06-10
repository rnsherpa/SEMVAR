#!/bin/sh

# Add a user with the same UID and GID as the host user
if [ -z "$HOST_UID" ]; then
    echo "HOST_UID is not set. Please pass your user UID."
    exit 1
fi

if [ -z "$HOST_GID" ]; then
    echo "HOST_GID is not set. Please pass your user GID."
    exit 1
fi

# Create a group with the HOST_GID
if ! getent group semgroup; then
    addgroup --gid $HOST_GID semgroup
fi

# Create a user with the HOST_UID
if ! id -u semuser >/dev/null 2>&1; then
    adduser --disabled-password --gecos '' --uid $HOST_UID --gid $HOST_GID semuser
fi

chown -R semuser:semgroup /SEMVAR

exec su semuser -c "python SEMVAR.py $@"