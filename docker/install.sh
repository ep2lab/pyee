#!/bin/bash
# This scripts makes your user and group IDs match inside the container so you can later modify any file created by the container

docker run -it -d --name setIDs -v $PYEE/docker:/root ep2lab/pyee:latest     # Run a detached interactive container
docker exec setIDs /bin/bash set_ids.sh `id -u` `id -g`                      # Execute the script creating the new user and group
docker commit setIDs ep2lab/pyee:local                                       # Commit to a new local image (this is run by pyee scripts)
docker stop setIDs                                                           # Stop and commit the container       
docker container rm setIDs