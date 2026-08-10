#!/bin/bash
# Must use bash for now due to substitution

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

JellyFish="$RUFUS_ROOT/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish"
Jhash="$1" # File to count within
HashList="$2" # List of kmers
MinCov="$3"
MaxCov="$4"

# This searches the HashList for specific kmers listed in the Jhash argument, and counts them
tmp=$(mktemp)
cat "$HashList" > "$tmp"
timeout 1h "$JellyFish" query -s <(awk '{print ">"$1"\n"$1}' "$tmp") "$Jhash" | awk -v var="$MinCov" ' $2 >= var ' | awk -v var="$MaxCov" ' $2 <= var '
rm -f "$tmp"
# timeout 1h "$JellyFish" query -s <(cat "$HashList" | awk '{print ">"$1"\n"$1}') "$Jhash" | awk -v var="$MinCov" ' $2 >= var ' | awk -v var="$MaxCov" ' $2 <= var '
#cat $HashList | awk '{print $1}' | $JellyFish query -i $Jhash | awk -v var=$MinCov ' $2 >= var ' | awk -v var=$MaxCov ' $2 <= var '
