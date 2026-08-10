#!/bin/bash

sudo mount /dev/nvme1n1 /mnt/data
sudo chown -R ubuntu:ubuntu /mnt/data

# if need to start from zero with new volume

#sudo file -s /dev/nvme1n1
#sudo mkfs -t ext4 /dev/nvme1n1
#sudo mkdir -p /mnt/data
#sudo mount /dev/nvme1n1 /mnt/data
#df -h
#sudo blkid /dev/nvme1n1
#echo "UUID=<your-uuid-here>  /mnt/data  ext4  defaults,nofail  0  2" | sudo tee -a /etc/fstab
#sudo chown -R ubuntu:ubuntu /mnt/data
